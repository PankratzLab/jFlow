// -Xmx4g
// -Xmx1024M
// found b132 in ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/
package bioinformatics;

import filesys.SerialHash;
import filesys.SnpMarkerSet;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import common.Aliases;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ProgressMonitor;
import common.ext;

public class ParseSNPlocations {
	public static final String DEFAULT_B36_SOURCE_FILENAME = "b130_SNPChrPosOnRef_36_3.bcp.gz";
	public static final String DEFAULT_B36_DB_FILENAME = "b130_36_3.ser";

	public static final String DEFAULT_B37_SOURCE_FILENAME = "b138_SNPChrPosOnRef.bcp.gz";
	public static final String DEFAULT_B37_DB_FILENAME = "b138_37_3.ser";

	public static final String DEFAULT_MERGE_SOURCE_FILENAME = "RsMergeArch.bcp.gz";
	public static final String DEFAULT_MERGE_FILENAME = "RsMerge.ser";
	
	public static final String DEFAULT_B37_VCF_FILENAME = "All_20150605_rsDominant.vcf.gz";
	public static final String DEFAULT_MERGE_VCF_FILENAME = "RsMerge.vcf.gz";
	public static final String DEFAULT_UNMAPPED_VCF_FILENAME = "All_20150605_papu_rsDominant.vcf.gz";
	
	public static final int DEFAULT_BUILD = 37;
	
	public static final int OFFSET = 1;
	public static final int MULTIPLE_POSITIONS = -2;
	public static final int UNMAPPED = -9;

//	public static void lowMemParse(String snpListFile, boolean useExistingPositions, Logger log) {
//		lowMemParse(snpListFile, ParseSNPlocations.DEFAULT_B37_DB, ParseSNPlocations.DEFAULT_MERGE, useExistingPositions, log);
//	}

	private static String[] TAG_SET = {"NSF","NSM","NSN","SYN","U3","U5","ASS","DSS","INT","R3","R5"};
	private static String[] TAG_SET_2 = {"PM","MUT"};
	
	public static int[][] parseSNPLocations(String snpListFile, String vcfFile, String unmappedVCF, String mergedVCF, Logger log, ProgressMonitor monitor) {
	    VCFFileReader vcfReader, unmappedVCFReader, mergedVCFReader;
        BufferedReader reader;
        
        ArrayList<int[]> resultList = new ArrayList<int[]>();
        
        String[] parts;
        String line = null;
        ArrayList<String> nonRS = new ArrayList<String>();
        ArrayList<String> rsNotFound = new ArrayList<String>();
        HashMap<String, Integer> indexMap = new HashMap<String, Integer>(); 
        int index = 0;
        
        String PROG_KEY = "PARSESNPS";
        int lineCnt = Files.countLines(snpListFile, 0);
        if (monitor != null) {
            monitor.beginDeterminateTask(PROG_KEY, "Processing " + lineCnt + " SNPs", lineCnt, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
        } else {
            System.out.println("Processing " + lineCnt + " SNPs");
        }
        
        vcfReader = new VCFFileReader(vcfFile, true);
        unmappedVCFReader = unmappedVCF == null ? null : new VCFFileReader(unmappedVCF, true);
        mergedVCFReader = mergedVCF == null ? null : new VCFFileReader(mergedVCF, true);
        try {
            reader = Files.getAppropriateReader(snpListFile);

            while ((line = reader.readLine()) != null) {
                parts = line.trim().split("[\\s]+");
                if (parts[0].equalsIgnoreCase("snp") || parts[0].toLowerCase().startsWith("marker")) {
                    continue;
                }
                index++;
                if (!parts[0].startsWith("rs")) {
                    resultList.add(new int[]{Integer.parseInt(parts.length > 1 ? parts[1] : "-1"),  Integer.parseInt(parts.length > 2 ? parts[2] : "-1")});
                    nonRS.add(parts[0]);
                    indexMap.put(parts[0], index);
                    // parse chr:pos:alleles markers later
                } else {
                    String rs = parts[0];
                    if (rs.indexOf(":") != -1) {
                        rs = rs.substring(0, rs.indexOf(":")); // remove alleles if they exist
                    }
                    int rsNumber = Integer.parseInt(rs.substring(2));// error checking on RS w/ allele
                    int chrom = rsNumber / (512 * 1024 * 1024) + 1;
                    
                    CloseableIterator<VariantContext> vcIter = vcfReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                    VariantContext markerVC = null;
                    while (vcIter.hasNext()) {
                        VariantContext vc = vcIter.next();
                        if (vc.getID().equals(rs)) {
                            markerVC = vc;
                            break;
                        }
                    }
                    vcIter.close();
                    vcIter = null;
                    if (markerVC == null) {
//                        System.err.println("Error - couldn't find {" + rs + "} in the regular database.  Checking merged and unmapped marker databases...");
                        if (unmappedVCFReader != null) {
                            CloseableIterator<VariantContext> vcIterUn = unmappedVCFReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                            while (vcIterUn.hasNext()) {
                                VariantContext vc = vcIterUn.next();
                                if (vc.getID().equals(rs)) {
                                    markerVC = vc;
                                    break;
                                }
                            }
                            vcIterUn.close();
                            vcIterUn = null;
                            if (markerVC == null) {
//                                System.err.println("Error - couldn't find {" + parts[0] + "} in the unmapped database.  Checking merged database now...");
                            }
                        }
                        if (markerVC == null && mergedVCFReader != null) {
                            String curr = null;
                            String next = rs;
                            HashSet<String> found = new HashSet<String>();
                            found.add(rs);
                            do {
                                curr = next;
                                next = null;
                                int rsCurr = Integer.parseInt(curr.startsWith("rs") ? curr.substring(2) : curr);
                                int chromCurr = rsCurr / (512 * 1024 * 1024) + 1;
                                CloseableIterator<VariantContext> vcIterUn = mergedVCFReader.query("chr" + chromCurr, rsCurr-2, rsCurr+2);
                                while (vcIterUn.hasNext()) {
                                    VariantContext vc = vcIterUn.next();
                                    if (vc.getID().equals(curr)) {
                                        if (!found.contains("rs" + vc.getAttribute("RSMRG").toString())) {
                                            next = "rs" + vc.getAttribute("RSMRG").toString();
                                            found.add("rs" + vc.getAttribute("RSMRG").toString());
                                        }
                                        break;
                                    }
                                }
                                vcIterUn.close();
                                vcIterUn = null;
                            } while (next != null);
                            
                            if (!curr.equals(rs)) {
                                rsNumber = Integer.parseInt(curr.substring(2));// error checking on RS w/ allele
                                chrom = rsNumber / (512 * 1024 * 1024) + 1;
                                vcIter = vcfReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                                markerVC = null;
                                while (vcIter.hasNext()) {
                                    VariantContext vc = vcIter.next();
                                    if (vc.getID().equals(curr)) {
                                        markerVC = vc;
                                        break;
                                    }
                                }
                                vcIter.close();
                                vcIter = null;
                            } else {
//                                System.err.println("Error - couldn't find {" + rs + "} in the merged database.");
                            }
                        }
                    }
                    if (markerVC != null) {
                        String attr = (String) markerVC.getAttribute("CHRPOS");
                        String[] pts = attr.split(":");
                        resultList.add(new int[]{Integer.parseInt(parts.length > 1 ? parts[1] : pts[0]), Integer.parseInt(parts.length > 2 ? parts[2] : pts[1])});
                    } else {
                        resultList.add(new int[]{Integer.parseInt(parts.length > 1 ? parts[1] : "-1"), Integer.parseInt(parts.length > 2 ? parts[2] : "-1")});
                        rsNotFound.add(rs);
                        indexMap.put(rs, index);
                    }
                }
                if (monitor != null) {
                    monitor.updateTask(PROG_KEY);
                }
            }
            reader.close();
        } catch (NumberFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        vcfReader.close();
        if (unmappedVCFReader != null) {
            unmappedVCFReader.close();
        }
        if (mergedVCFReader != null) {
            mergedVCFReader.close();
        }
        if (rsNotFound.size() > 0) {
            Files.writeArrayList(rsNotFound, ext.rootOf(snpListFile, false)+"_missing.txt");
        }
        if (monitor != null) {
            monitor.endTask(PROG_KEY);
        }
        
        return resultList.toArray(new int[resultList.size()][]);
	}
	
	public static void parseSNPlocations(String snpListFile, String vcfFile, String unmappedVCF, String mergedVCF, Logger log, ProgressMonitor monitor) {
        VCFFileReader vcfReader, unmappedVCFReader, mergedVCFReader;
        BufferedReader reader;
        PrintWriter writer;
        
        String[] parts;
        String line = null;
        ArrayList<String> nonRS = new ArrayList<String>();
        ArrayList<String> rsNotFound = new ArrayList<String>();
        HashMap<String, Integer> indexMap = new HashMap<String, Integer>(); 
        int index = 0;
        
        String PROG_KEY = "PARSESNPS";
        int lineCnt = Files.countLines(snpListFile, 0);
        if (monitor != null) {
            monitor.beginDeterminateTask(PROG_KEY, "Processing " + lineCnt + " SNPs", lineCnt, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
        } else {
            System.out.println("Processing " + lineCnt + " SNPs");
        }
        
        vcfReader = new VCFFileReader(vcfFile, true);
        unmappedVCFReader = unmappedVCF == null ? null : new VCFFileReader(unmappedVCF, true);
        mergedVCFReader = mergedVCF == null ? null : new VCFFileReader(mergedVCF, true);
        try {
            reader = Files.getAppropriateReader(snpListFile);
            writer = Files.getAppropriateWriter(ext.rootOf(snpListFile, false)+"_positions.xln");
            writer.println("SNP\tChr\tPosition\tRef\tAlt\tFunc\tPM/MUT\tGENENAME\tCAF");

            while ((line = reader.readLine()) != null) {
                parts = line.trim().split("[\\s]+");
                if (parts[0].equalsIgnoreCase("snp") || parts[0].toLowerCase().startsWith("marker")) {
                    continue;
                }
                index++;
                if (!parts[0].startsWith("rs")) {
                    writer.println(parts[0] + "\t" + (parts.length > 1 ? parts[1] : ".") + "\t" + (parts.length > 2 ? parts[2] : ".") + "\t" + (parts.length > 3 ? parts[3] : ".") + "\t" + (parts.length > 4 ? parts[4] : ".") + "\t" + "." + "\t" + "." + "\t" + "." + "\t" + ".");
                    nonRS.add(parts[0]);
                    indexMap.put(parts[0], index);
                    // parse chr:pos:alleles markers later
                } else {
                    String rs = parts[0];
                    if (rs.indexOf(":") != -1) {
                        rs = rs.substring(0, rs.indexOf(":")); // remove alleles if they exist
                    }
                    int rsNumber = Integer.parseInt(rs.substring(2));// error checking on RS w/ allele
                    int chrom = rsNumber / (512 * 1024 * 1024) + 1;
                    
                    CloseableIterator<VariantContext> vcIter = vcfReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                    VariantContext markerVC = null;
                    while (vcIter.hasNext()) {
                        VariantContext vc = vcIter.next();
                        if (vc.getID().equals(rs)) {
                            markerVC = vc;
                            break;
                        }
                    }
                    vcIter.close();
                    vcIter = null;
                    if (markerVC == null) {
//                        System.err.println("Error - couldn't find {" + rs + "} in the regular database.  Checking merged and unmapped marker databases...");
                        if (unmappedVCFReader != null) {
                            CloseableIterator<VariantContext> vcIterUn = unmappedVCFReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                            while (vcIterUn.hasNext()) {
                                VariantContext vc = vcIterUn.next();
                                if (vc.getID().equals(rs)) {
                                    markerVC = vc;
                                    break;
                                }
                            }
                            vcIterUn.close();
                            vcIterUn = null;
                            if (markerVC == null) {
//                                System.err.println("Error - couldn't find {" + parts[0] + "} in the unmapped database.  Checking merged database now...");
                            }
                        }
                        if (markerVC == null && mergedVCFReader != null) {
                            String curr = null;
                            String next = rs;
                            HashSet<String> found = new HashSet<String>();
                            found.add(rs);
                            do {
                                curr = next;
                                next = null;
                                int rsCurr = Integer.parseInt(curr.startsWith("rs") ? curr.substring(2) : curr);
                                int chromCurr = rsCurr / (512 * 1024 * 1024) + 1;
                                CloseableIterator<VariantContext> vcIterUn = mergedVCFReader.query("chr" + chromCurr, rsCurr-2, rsCurr+2);
                                while (vcIterUn.hasNext()) {
                                    VariantContext vc = vcIterUn.next();
                                    if (vc.getID().equals(curr)) {
                                        if (!found.contains("rs" + vc.getAttribute("RSMRG").toString())) {
                                            next = "rs" + vc.getAttribute("RSMRG").toString();
                                            found.add("rs" + vc.getAttribute("RSMRG").toString());
                                        }
                                        break;
                                    }
                                }
                                vcIterUn.close();
                                vcIterUn = null;
                            } while (next != null);
                            
                            if (!curr.equals(rs)) {
                                rsNumber = Integer.parseInt(curr.substring(2));// error checking on RS w/ allele
                                chrom = rsNumber / (512 * 1024 * 1024) + 1;
                                vcIter = vcfReader.query("chr" + chrom, rsNumber-2, rsNumber+2);
                                markerVC = null;
                                while (vcIter.hasNext()) {
                                    VariantContext vc = vcIter.next();
                                    if (vc.getID().equals(curr)) {
                                        markerVC = vc;
                                        break;
                                    }
                                }
                                vcIter.close();
                                vcIter = null;
                            } else {
//                                System.err.println("Error - couldn't find {" + rs + "} in the merged database.");
                            }
                        }
                    }
                    if (markerVC != null) {
                        String attr = (String) markerVC.getAttribute("CHRPOS");
                        String[] pts = attr.split(":");
                        StringBuilder newLine = new StringBuilder();
                        newLine.append(rs).append("\t")
                                .append(parts.length > 1 ? parts[1] : pts[0]).append("\t")
                                .append(parts.length > 2 ? parts[2] : pts[1]).append("\t");
                        String ref = markerVC.getReference().toString();
                        if (ref.endsWith("*")) {
                            ref = ref.substring(0, ref.length() - 1);
                        }
                        newLine.append(ref).append("\t") 
                                .append(markerVC.getAltAlleleWithHighestAlleleCount()).append("\t");
                        boolean found = false;
                        for (String funcKey : TAG_SET) {
                            if (markerVC.getAttributes().keySet().contains(funcKey)) {
                                if (found) {
                                    newLine.append(";");
                                } else {
                                    found = true;
                                }
                                newLine.append(funcKey);
                            }
                        }
                        if (!found) {
                            newLine.append(".");
                        }
                        newLine.append("\t");
                        found = false;
                        for (String funcKey : TAG_SET_2) {
                            if (markerVC.getAttributes().keySet().contains(funcKey)) {
                                if (found) {
                                    newLine.append(";");
                                } else {
                                    found = true;
                                }
                                newLine.append(funcKey);
                            }
                        }
                        if (!found) {
                            newLine.append(".");
                        }
                        newLine.append("\t");
                        java.util.Map<String, Object> attrs = markerVC.getAttributes();
                        if (attrs.containsKey("GENEINFO")) {
                            newLine.append(attrs.get("GENEINFO")).append("\t");
                        } else {
                            newLine.append(".\t");
                        }
                        ArrayList<?> caf = (ArrayList<?>) attrs.get("CAF");
                        if (caf == null) {
                            newLine.append(".");
                        } else {
                            for (int i = 0; i < caf.size(); i++) {
                                newLine.append(caf.get(i));
                                if (i < caf.size() - 1) {
                                    newLine.append(";");
                                }
                            }
                        }
                        writer.println(newLine.toString());
                    } else {
                        writer.println(rs + "\t" + (parts.length > 1 ? parts[1] : ".") + "\t" + (parts.length > 2 ? parts[2] : ".") + "\t" + (parts.length > 3 ? parts[3] : ".") + "\t" + (parts.length > 4 ? parts[4] : ".") + "\t" + "." + "\t" + "." + "\t" + "." + "\t" + ".");
                        rsNotFound.add(rs);
                        indexMap.put(rs, index);
                    }
                }
                if (monitor != null) {
                    monitor.updateTask(PROG_KEY);
                }
            }
            writer.flush();
            writer.close();
            reader.close();
        } catch (NumberFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        vcfReader.close();
        if (unmappedVCFReader != null) {
            unmappedVCFReader.close();
        }
        if (mergedVCFReader != null) {
            mergedVCFReader.close();
        }
        if (rsNotFound.size() > 0) {
            Files.writeArrayList(rsNotFound, ext.rootOf(snpListFile, false)+"_missing.txt");
        }
        if (monitor != null) {
            monitor.endTask(PROG_KEY);
        }
    }
	
	public static void lowMemParse(String snpListFile, String db, String mergeDB, boolean useExistingPositions, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		SnpMarkerSet dbMarkerSet;
		int[] dbRSnumbers, dbPositions;
		byte[] dbChrs;
		int index;
		Hashtable<Integer,Integer> mergeHash;
		Integer trav, next;
		boolean first;
		byte chr;
		int position;

		log.report(ext.getTime());

		dbMarkerSet = null;
		dbRSnumbers = null;
		dbPositions = null;
		dbChrs = null;
		mergeHash = null;
		try {
	        reader = new BufferedReader(new FileReader(snpListFile));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(snpListFile, false)+"_positions.xln"));
			writer.println("SNP\tChr\tPosition");
			first = true;
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	chr = -2;
	        	position = -2;
	        	if (first) {
	        		if (line[0].equalsIgnoreCase("snp") || line[0].toLowerCase().startsWith("marker")) {
	    	        	line = reader.readLine().trim().split("[\\s]+");
	        		}
	        		first = false;
	        	}
	        	if (line.length > 2 && !line[1].equals(".") && !line[2].equals(".")) {
					writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
	        	} else {
	        		if (dbMarkerSet == null) {
	        			log.report("Loading database...");
	        			log.report("(if file is not found and it's searching for the wrong file, then force a recompile of MapSNPsAndGenes so that it gets the new constant from ParseSNPlocations)");

	        			dbMarkerSet = SnpMarkerSet.load(db, false, log);
	        			dbRSnumbers = dbMarkerSet.getRSnumbers();
	        			dbPositions = dbMarkerSet.getPositions();
	        			dbChrs = dbMarkerSet.getChrs();

	        			log.report("Searching database...");
	        			
	        		}
	        		if (line[0].startsWith("rs")) {
	        		    String rs = line[0];
	        		    if (rs.indexOf(":") != -1) {
	        		        rs = rs.substring(0, rs.indexOf(":")); // remove alleles if they exist
	        		    }
						index = Array.binarySearch(dbRSnumbers, Integer.parseInt(rs.substring(2)), true); // TODO parse RS w/ alleles
						if (index == -1) {
							if (mergeDB != null) {
								if (mergeHash == null) {
									log.report("Could not find a marker, loading merge database...", false, true);
									if (new File(mergeDB).exists()) {
										mergeHash = SerialHash.loadSerializedIntHash(mergeDB);
										log.report("done");
									} else {
										log.report("failed; "+mergeDB+" not found");
										mergeHash = new Hashtable<Integer,Integer>();
									}
								}
								trav = Integer.parseInt(rs.substring(2));
								next = trav;
								while (next != null) {
									trav = next;
									next = mergeHash.get(trav);
								}
								if (rs.equals("rs"+trav)) {
									log.reportError("\n\n****ERROR**** failed to find "+rs+" in any NCBI database ****ERROR****\n\n");
								} else if (trav != null) {
									log.reportError("FYI - "+rs+" has merged with rs"+trav);
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
									
									if (index == -1) {
										log.reportError("Error - could not find rs"+trav+" in "+db+"; must have been added after this db was released");
										chr = (byte)0;
										position = 0;
									} else {
										chr = dbChrs[index];
										position = dbPositions[index]+ParseSNPlocations.OFFSET;
										index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
									}
								} else {
									log.reportError("Error - tried and failed to find a merge record for "+rs);
								}
							}
							if (index == -1) {
								log.reportError("Error - could not find "+rs+" in "+db);
								chr = (byte)0;
								position = 0;
							}
						} else {
							if (dbChrs[index] == ParseSNPlocations.UNMAPPED) {
								log.reportError("Error - marker "+rs+" is not mapped to a chromosome");
								chr = 0;
								position = -1;
							} else if (dbPositions[index] == ParseSNPlocations.MULTIPLE_POSITIONS) {
								log.reportError("Warning - marker "+rs+" likely has multiple positions on chromosome "+dbChrs[index]);
								chr = dbChrs[index];
								position = -1;
							} else {
								chr = dbChrs[index];
								position = dbPositions[index]+ParseSNPlocations.OFFSET;
							}
						}
						writer.println(rs+"\t"+chr+"\t"+position);
	        		} else if (line[0].startsWith("chr") || line[0].indexOf(":") != -1) {
	        		    String[] pts = line[0].split(":");
	        		    int chrTemp = -2;
	        		    int posTemp = -2;
	        		    if (pts.length >= 2) {
    	        		    try {
    	        		    	pts[0] = ext.replaceAllWith(pts[0], new String[][] {{"X", "23"}, {"Y", "24"}, {"XY", "25"}, {"M", "26"}});
    	        		        chrTemp = Integer.parseInt(pts[0].startsWith("chr") ? pts[0].substring(3) : pts[0]);
    	        		    } catch (NumberFormatException e) {}
    	        		    try {
    	        		        posTemp = Integer.parseInt(pts[1]);
    	        		    } catch (NumberFormatException e) {}
	        		    }
	        		    if (chrTemp == -2 || posTemp == -2) {
	                        log.reportError("Error - can't look up a SNP without an rs number or chromosome/position name ("+line[0]+")");
	                        writer.println(line[0]+"\t0\t0");
	        		    } else {
	                        writer.println(line[0]+"\t" + chrTemp + "\t" + posTemp);
	        		    }
					} else {
						log.reportError("Error - can't look up a SNP without an rs number or chromosome/position name ("+line[0]+")");
						writer.println(line[0]+"\t0\t0");
					}
	        	}
	        }
	        reader.close();
	        writer.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+snpListFile+"\" not found in current directory");
//	        System.exit(1);
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+snpListFile+"\"");
//	        System.exit(2);
        }

		log.report("Done with the positions! "+ext.getTime());
	}

	public static void hashParse(String snpListFile, String db, String mergeDB, boolean useExistingPositions, Logger log) {
//		boolean reminderToCode = true;
		// TODO Auto-generated catch block
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		SnpMarkerSet dbMarkerSet;
		int[] dbRSnumbers, dbPositions;
		byte[] dbChrs;
		int index;
		Hashtable<Integer,Integer> mergeHash;
		Integer trav, next;
		boolean first;
		byte chr;
		int position;

		log.report(ext.getTime());

		dbMarkerSet = null;
		dbRSnumbers = null;
		dbPositions = null;
		dbChrs = null;
		mergeHash = null;
		try {
	        reader = new BufferedReader(new FileReader(snpListFile));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(snpListFile, false)+"_positions.xln"));
			writer.println("SNP\tChr\tPosition");
			first = true;
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	chr = -2;
	        	position = -2;
	        	if (first) {
	        		if (line[0].equalsIgnoreCase("snp") || line[0].toLowerCase().startsWith("marker")) {
	    	        	line = reader.readLine().trim().split("[\\s]+");
	        		}
	        		first = false;
	        	}
	        	if (line.length > 2 && !line[1].equals(".") && !line[2].equals(".")) {
					writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
	        	} else {
	        		if (dbMarkerSet == null) {
	        			log.report("Loading database...");

	        			dbMarkerSet = SnpMarkerSet.load(db, false);
	        			dbRSnumbers = dbMarkerSet.getRSnumbers();
	        			dbPositions = dbMarkerSet.getPositions();
	        			dbChrs = dbMarkerSet.getChrs();

	        			log.report("Searching database...");
	        			
	        		}
	        		if (line[0].startsWith("rs")) {
						index = Array.binarySearch(dbRSnumbers, Integer.parseInt(line[0].substring(2)), true);
						if (index == -1) {
							if (mergeDB != null) {
								if (mergeHash == null) {
									log.report("Could not find a marker, loading merge database...", false, true);
									if (new File(mergeDB).exists()) {
										mergeHash = SerialHash.loadSerializedIntHash(mergeDB);
										log.report("done");
									} else {
										log.report("failed; "+mergeDB+" not found");
										mergeHash = new Hashtable<Integer,Integer>();
									}
								}
								trav = Integer.parseInt(line[0].substring(2));
								next = trav;
								while (next != null) {
									trav = next;
									next = mergeHash.get(trav);
								}
								if (line[0].equals("rs"+trav)) {
									log.reportError("\n\n****ERROR**** failed to find "+line[0]+" in any NCBI database ****ERROR****\n\n");
								} else if (trav != null) {
									log.reportError("FYI - "+line[0]+" has merged with rs"+trav);
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
									chr = dbChrs[index];
									position = dbPositions[index]+ParseSNPlocations.OFFSET;
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
								} else {
									log.reportError("Error - tried and failed to find a merge record for "+line[0]);
								}
							}
							if (index == -1) {
								log.reportError("Error - could not find "+line[0]+" in "+db);
								chr = (byte)0;
								position = 0;
							}
						} else {
							if (dbChrs[index] == ParseSNPlocations.UNMAPPED) {
								log.reportError("Error - marker "+line[0]+" is not mapped to a chromosome");
								chr = 0;
								position = -1;
							} else if (dbPositions[index] == ParseSNPlocations.MULTIPLE_POSITIONS) {
								log.reportError("Warning - marker "+line[0]+" likely has multiple positions on chromosome "+dbChrs[index]);
								chr = dbChrs[index];
								position = -1;
							} else {
								chr = dbChrs[index];
								position = dbPositions[index]+ParseSNPlocations.OFFSET;
							}
						}
						writer.println(line[0]+"\t"+chr+"\t"+position);
					} else {
						log.reportError("Error - can't look up a SNP without an rs number ("+line[0]+")");
						writer.println(line[0]+"\t0\t0");
					}
	        	}
	        }
	        reader.close();
	        writer.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+snpListFile+"\" not found in current directory");
//	        System.exit(1);
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+snpListFile+"\"");
//	        System.exit(2);
        }

		log.report("Done with the positions! "+ext.getTime());
	}
	
	public static void createFromSource(String source, String db) {
		BufferedReader reader;
		String[] line;
		int count, countPARys;
		String[] rsNumbers;
		int[] positions;
		byte[] chrs;
		String temp;

		try {
			System.out.println("Reading in "+source);
//			reader = new BufferedReader(new FileReader(source));
			reader = Files.getAppropriateReader(source);
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
//				if (line.length>2) {
				if (line.length>1) {
					count++;
				}
			}
			reader.close();
			System.out.println("Found "+count+" markers");

//			reader = new BufferedReader(new FileReader(source));
			reader = Files.getAppropriateReader(source);
			rsNumbers = new String[count];
			positions = new int[count];
			chrs = new byte[count];
			count = 0;
			countPARys = 0;
			while (reader.ready()) {
				temp = reader.readLine();
				try {
					line = temp.trim().split("[\\s]+");
//					if (line.length>2) {
					if (line.length>1) {
						rsNumbers[count] = "rs"+line[0];
						chrs[count] = Positions.chromosomeNumber(line[1]);
						if (line[1].equals("PAR") && line[2].equals("y")) {
							countPARys++;
						} else if (line.length > 2) { // new
							try {
								positions[count] = Integer.parseInt(line[2]);
							} catch (NumberFormatException nfe) {
								System.err.println("Error - parsing line #"+(count+1)+": "+Array.toStr(line));
							}
						} else {
							positions[count] = MULTIPLE_POSITIONS;
						}
						count++;
					} else { // new
						chrs[count] = 0;
						positions[count] = Integer.parseInt(line[2]);
						count++;
					}
				} catch (Exception e) {
					System.err.println("Error parsing "+source+" at the following line:");
					System.err.println(temp);
					e.printStackTrace();
					reader.close();
					return;
				}
				
			}
			reader.close();
			if (countPARys > 0) {
				System.err.println("Warning - there were "+countPARys+" instances where the chromosome was listed as PAR and the position was nothing but 'y'");
			}
			// already in order, d'uh!
//			System.out.println("Sorting...");
//			keys = Sort.quicksort(rsNumbers);
			System.out.println("Serializing...");
			new SnpMarkerSet(rsNumbers, chrs, positions, null, null, false, true).serialize(db);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+source+"\" not found in current directory");
//			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+source+"\"");
//			System.exit(2);
		}
	}

	public static void createMergeDBfromSource(String source, String db) {
		BufferedReader reader;
		String[] line;
		int count;
		Hashtable<Integer,Integer> hash;

		count = 0;
		try {
			System.out.println(ext.getTime()+"\tReading in "+source);
			reader = new BufferedReader(new FileReader(source));
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();
			System.out.println(ext.getTime()+"\tFound "+count+" merged markers");
			
			hash = new Hashtable<Integer,Integer>(count);

			count = 0;
			reader = Files.getAppropriateReader(source);
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\t", -1);
//				hash.put(line[0], line[6]); // not complete
				hash.put(Integer.parseInt(line[0]), Integer.parseInt(line[1]));
				count++;
			}
			reader.close();
			System.out.println(ext.getTime()+"\tSerializing...");
			SerialHash.createSerializedIntHash(db, hash);
			System.out.println(ext.getTime()+"\tDone!");
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+source+"\" not found in current directory");
//			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+source+"\"");
//			System.exit(2);
		} catch (OutOfMemoryError oome) {
			System.err.println("Error - Ran out of memory at line "+(count+1));
			oome.printStackTrace();
		}
	}	
	
	
	
	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String source = "";
		String mergeSource = "";
		Logger log = new Logger();
		String dir = Files.firstDirectoryThatExists(Aliases.REFERENCE_FOLDERS, false, false, log);
		String db = Aliases.getPathToFileInReferenceDirectory(DEFAULT_B37_DB_FILENAME, false, log);
		String merge = Aliases.getPathToFileInReferenceDirectory(DEFAULT_MERGE_FILENAME, false, log);
		String vcf = null;
		String mergedvcf = null;
		String unmappedvcf = null;

		// uncomment one of these to compile
//		source = DEFAULT_B37_SOURCE;
//		db = DEFAULT_B37_DB;

//		source = DEFAULT_B36_SOURCE;
//		db = DEFAULT_B36_DB;
		
//		dir = "";
//		mergeSource = "C:/Users/cole0482/Downloads/RsMergeArch.bcp.gz";
//		db = "D:/RSMerge.ser";
//		 vcf=N:/statgen/NCBI/All_20150605_rsDominant.vcf.gz unmappedvcf=N:/statgen/NCBI/All_20150605_papu_rsDominant.vcf.gz mergedvcf=N:/statgen/NCBI/RsMerge.vcf.gz
		
		String filename = "list.txt";
		boolean plinkFormat = false;

		String usage = "\n"+
		"bioinformatics.ParseSNPlocations requires 0-1 arguments\n"+
		"   (1) file from which to create snp db (i.e. source="+DEFAULT_B37_SOURCE_FILENAME+" (not the default))\n"+
		"  OR:\n"+
		"   (1) file from which to create merge db (i.e. source="+DEFAULT_MERGE_SOURCE_FILENAME+" (not the default))\n"+
		"  OR:\n"+
		"   (1) snpDB (i.e. db="+db+" (default))\n"+
		"   (2) merge DB (i.e. merge="+merge+" (default))\n"+
		"   (3) dir (i.e. dir="+dir+" (default))\n"+
		"   (4) filename (i.e. file="+filename+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("db=")) {
				db = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vcf=")) {
			    vcf = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("mergedvcf=")) {
			    mergedvcf = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("unmappedvcf=")) {
			    unmappedvcf = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("merge=")) {
				merge = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (!source.equals("")) {
				createFromSource(dir+source, dir+db);
			} else if (!mergeSource.equals("")) {
				createMergeDBfromSource(dir+mergeSource, dir+merge);
			} else {
				if (vcf != null) {
				    long t = System.currentTimeMillis();
				    log = new Logger();
				    parseSNPlocations(dir+filename, vcf, unmappedvcf, mergedvcf, log, new ProgressMonitor(null, log));
				    System.out.println("Took " + ext.getTimeElapsed(t));
				} else {
				    long t = System.currentTimeMillis();
				    SnpMarkerSet map = new SnpMarkerSet(dir+filename, plinkFormat?SnpMarkerSet.PLINK_MAP_FORMAT:SnpMarkerSet.NAMES_ONLY, true, new Logger());
				    map.parseSNPlocations(db, merge, new Logger());
				    map.writeToFile(dir+ext.rootOf(filename)+"_newPositions.out", SnpMarkerSet.GENERIC_FORMAT);
				    System.out.println("Took " + ext.getTimeElapsed(t));
				}
				
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
