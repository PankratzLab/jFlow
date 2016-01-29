// -Xms6G -Xmx6G
// expand to dynamically load/save a certain chunk of markers at a time  
package filesys;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinformatics.Sequence;
import common.*;
import stats.*;
import gwas.*;

public class DosageData implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final int MACH_MLDOSE_FORMAT = 0;
	public static final int GEN_FORMAT = 1;
	public static final int GWAF_FORMAT = 2;
	public static final int PLINK_FORMAT = 3;
	public static final int MACH_MLPROB_FORMAT = 4;
	public static final int MINIMAC_DOSE_FORMAT = 5;
	public static final int IMPUTE2_DOSE_FORMAT = 6;
	public static final int DATABASE_DOSE_FORMAT = 7;
	public static final int BEAGLE_DOSE_FORMAT = 8; // not currently included in determineType

	public static final int MACH_ID_TYPE = 0;
	public static final int SEPARATE_FILE_ID_TYPE = 1;
	public static final int IID_TYPE = 2;
	public static final int FID_IID_TYPE = 3;
	
	public static final int INDIVIDUAL_DOMINANT_FORMAT = 0;
	public static final int MARKER_DOMINANT_FORMAT = 1;
	
	public static final int CHR_INFO_IN_FILENAME = -2;
	
	public static final String CHR_REGEX = ".*?chr(\\d\\d?).*?";
	
	public static final String[][] HEADS = {null,				null,	{"id"}, 	{"SNP", "A1", "A2"},	null, 				null,				{"FID", "IID"}};
	public static final String[][] LEADS = {{null, "MLDOSE"},	null,	null, 		null, 					{null, "MLPROB"},	{null, "DOSE"},		null};
	public static final String[] DELIMITERS = {"\t", ",", " "};

	/** 0,       1,                                                  2,                3,                                      4,          5,                6,        7,        8,         9,         10,        11,                  12,                   13                   */
	/** id_type, column index where dosage/probability values begin, dominance format, number of columns summarizing the data, header row, marker/IID index, A1 index, A2 index, chr index, pos index, delimiter, index for head/lead, min number of digits, max number of digits */

	/** 									   0                      1  2                           3  4  5   6   7   8   9 10 11 12 13 */
	public static final int[][] PARAMETERS = {{MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1, -1, -1, -1, 0, 0, 3, 3}, // .mldose (MACH)
											  {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT,     3, 0, 1,  3,  4,  0,  2, 2, 1, 0, 3}, // .gen
											  {IID_TYPE,              1, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0, -1, -1, -1, -1, 1, 2, 0, 3}, // .fhsR (GWAF) 
											  {FID_IID_TYPE,          3, MARKER_DOMINANT_FORMAT,     2, 1, 0,  1,  2, -1, -1, 0, 3, 3, 3}, // .dosage (PLINK)
											  {MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 2, 0, 0, -1, -1, -1, -1, 0, 4, 3, 3}, // .mlprob (MACH)
											  {MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1, -1, -1, -1, 0, 5, 3, 3}, // .dose (MINIMAC)
											  {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT,     3, 0, 1,  3,  4, CHR_INFO_IN_FILENAME,  2, 0, 1, 3, 3}, // .impute2
											  {FID_IID_TYPE,          3, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0,  3,  4, -1,  2, 0, 6, 3, 3}, // .db.xln
											  {SEPARATE_FILE_ID_TYPE, 3, MARKER_DOMINANT_FORMAT,     1, 1, 0,  1,  2, -1, -1, 2, 1, 4, 4}, // .dose (BEAGLE)
	};

	private SnpMarkerSet markerSet;
	private String[][] ids;
	private float[][] dosageValues;
	private float[][][] genotypeProbabilities;
	private char[][] alleles;
	private byte[] chrs;
	private int[] positions;
	
	public static DosageData combine(DosageData dd1, DosageData dd2) {
	    byte missingChr = 0;
	    int missingPos = 0;
	    char[] missingAlleles = null;
	    float missingDosage = 0;
	    float missingGeno = 0;
	    
	    
	    String[][] dd1Ids = dd1.ids;
	    String[][] dd2Ids = dd2.ids;
	    LinkedHashSet<String> idSet = new LinkedHashSet<String>(); // use to ensure uniqueness and order
	    for (String[] id : dd1Ids) {
	        idSet.add(id[0] + "\t" + id[1]);
	    }
	    HashSet<Integer> duplicatedSampleIndices = new HashSet<Integer>();
	    HashMap<String, Integer> duplicatedSamplesAndIndices = new HashMap<String, Integer>();
	    for (int i = 0; i < dd2Ids.length; i++) {
	        boolean alreadyPresent = !idSet.add(dd2Ids[i][0] + "\t" + dd2Ids[i][1]);
	        if (alreadyPresent) {
	            duplicatedSampleIndices.add(i);
	            duplicatedSamplesAndIndices.put(dd2Ids[i][0] + "\t" + dd2Ids[i][1], i);
	            // TODO log? mark positions?
	        }
	    }
	    if (duplicatedSampleIndices.size() > 0) {
	        System.out.println(duplicatedSampleIndices.size() + " duplicate sample IDs found, out of " + dd2Ids.length + " samples present.");
	    }
	    
	    String[] dd1Mkrs = dd1.markerSet.getMarkerNames();
	    String[] dd2Mkrs = dd2.markerSet.getMarkerNames();
	    LinkedHashSet<String> markers = new LinkedHashSet<String>();
	    for (String mkr : dd1Mkrs) {
	        markers.add(mkr);
	    }
	    HashSet<Integer> duplicatedMarkerIndices = new HashSet<Integer>();
	    HashMap<String, Integer> duplicatedMarkersAndIndices = new HashMap<String, Integer>();
	    HashMap<String, Integer> dd2markersAndIndices = new HashMap<String, Integer>();
	    for (int i = 0; i < dd2Mkrs.length; i++) {
	        dd2markersAndIndices.put(dd2Mkrs[i], i);
	        boolean alreadyPresentMkr = !markers.add(dd2Mkrs[i]);
	        if (alreadyPresentMkr) {
	            duplicatedMarkerIndices.add(i);
	            duplicatedMarkersAndIndices.put(dd2Mkrs[i], i);
	            if (duplicatedSampleIndices.size() > 0) {
	                System.out.println("Duplicate marker: " + dd2Mkrs[i]);
	                System.err.println("Error - cannot combine data sets with the same marker and sample names.  Yet.");
	                // TODO 
	                // combining values when marker and sample are the same?
	                System.exit(1);
	            }
	        }
	    }
	    
	    DosageData ddNew = new DosageData();
	    ddNew.ids = new String[idSet.size()][];
	    Iterator<String> iter = idSet.iterator();
	    int ind = 0;
	    while (iter.hasNext()) {
	        ddNew.ids[ind++] = iter.next().split("\t");
	    }
	    idSet = null; // can now refer to ddNew.ids
	    
	    ddNew.markerSet = SnpMarkerSet.merge(dd1.markerSet, dd2.markerSet);
	    
	    ddNew.alleles = new char[markers.size()][];
	    ddNew.chrs = new byte[markers.size()];
	    ddNew.positions = new int[markers.size()];
	    byte[] chrSrc = dd1.chrs == null ? dd1.markerSet.getChrs() : dd1.chrs;
	    char[][] alleleSrc = dd1.alleles == null ? dd1.markerSet.getAlleles() : dd1.alleles;
	    int[] posSrc = dd1.positions == null ? dd1.markerSet.getPositions() : dd1.positions;
	    int chrOffset = chrSrc == null ? dd1Mkrs.length : chrSrc.length;
	    for (int i = 0; i < chrSrc.length; i++) {
	        ddNew.chrs[i] = chrSrc == null ? missingChr : chrSrc[i];
	        ddNew.alleles[i] = alleleSrc == null ? missingAlleles : alleleSrc[i];
	        ddNew.positions[i] = posSrc == null ? missingPos : posSrc[i];
	    }
	    chrSrc = dd2.chrs == null ? dd2.markerSet.getChrs() : dd2.chrs;
	    alleleSrc = dd2.alleles == null ? dd2.markerSet.getAlleles() : dd2.alleles;
	    posSrc = dd2.positions == null ? dd2.markerSet.getPositions() : dd2.positions;
	    for (int i = 0; i < chrSrc.length; i++) {
	        if (!duplicatedMarkerIndices.contains(i)) {
	            ddNew.chrs[chrOffset + i] = chrSrc == null ? missingChr : chrSrc[i];
	            ddNew.alleles[chrOffset + i] = alleleSrc == null ? missingAlleles : alleleSrc[i];
	            ddNew.positions[chrOffset + i] = posSrc == null ? missingPos : posSrc[i];
	        }
	    }
	    
	    int dd1NumGeno = dd1.genotypeProbabilities == null ? 0 : dd1.genotypeProbabilities[0][0].length;
        int dd2NumGeno = dd2.genotypeProbabilities == null ? 0 : dd2.genotypeProbabilities[0][0].length;
	    
        if ((dd1NumGeno > 1 || dd2NumGeno > 1) && dd1NumGeno != dd2NumGeno) {
            System.out.println("dd1: " + dd1NumGeno + " | dd2: " + dd2NumGeno);
            System.err.println("Error - cannot combine data sets with different numbers of genotype probabilities.  Yet.");
            // TODO error, mismatched number of genotype probabilities
            System.exit(1);
        }
        int ddNewNumGeno = dd1NumGeno; 
        
        ddNew.genotypeProbabilities = ddNewNumGeno > 1 ? new float[markers.size()][ddNew.ids.length][ddNewNumGeno] : null;
        ddNew.dosageValues = ddNewNumGeno == 1 ? new float[markers.size()][ddNew.ids.length] : null;
        
        if (ddNewNumGeno > 1) {
            // combine genotypeProbs
            Iterator<String> markerIter = markers.iterator();
            int m = 0;
            int dd2IndM = -1;
            int dd1MarkerOffset = dd1Mkrs.length;
            while (markerIter.hasNext()) {
                String mkr = markerIter.next();
                dd2IndM = -1;
                if (duplicatedMarkersAndIndices.containsKey(mkr)) {
                    dd2IndM = duplicatedMarkersAndIndices.get(mkr);
                }
                
                for (int s = 0; s < dd1.ids.length; s++) {
                    if (m < dd1MarkerOffset) {
                        ddNew.genotypeProbabilities[m][s] = dd1.genotypeProbabilities[m][s];
                    } else {
                        if (duplicatedSamplesAndIndices.containsKey(dd1.ids[s][0] + "\t" + dd1.ids[s][1])) {
                            ddNew.genotypeProbabilities[m][s] = dd2.genotypeProbabilities[dd2markersAndIndices.get(mkr)][duplicatedSamplesAndIndices.get(dd1.ids[s][0] + "\t" + dd1.ids[s][1])];
                        } else {
                            ddNew.genotypeProbabilities[m][s] = Array.floatArray(ddNewNumGeno, missingGeno);
                        }
                    }
                }
                int newInd = dd1.ids.length;
                for (int s = 0; s < dd2.ids.length; s++) {
                    if (duplicatedSampleIndices.contains(s)) {
                        continue;
                    }
                    if (!dd2markersAndIndices.containsKey(mkr)) {
                        ddNew.genotypeProbabilities[m][newInd] = Array.floatArray(ddNewNumGeno, missingGeno);
                    } else {
                        ddNew.genotypeProbabilities[m][newInd] = dd2.genotypeProbabilities[dd2IndM == -1 ? dd2markersAndIndices.get(mkr) : dd2IndM][s];
                    }
                    newInd++;
                }
                
                m++;
            }
        } else if (ddNewNumGeno == 1) {

            Iterator<String> markerIter = markers.iterator();
            int m = 0;
            int dd2IndM = -1;
            int dd1MarkerOffset = dd1Mkrs.length;
            while (markerIter.hasNext()) {
                String mkr = markerIter.next();
                dd2IndM = -1;
                if (duplicatedMarkersAndIndices.containsKey(mkr)) {
                    dd2IndM = duplicatedMarkersAndIndices.get(mkr);
                }
                
                for (int s = 0; s < dd1.ids.length; s++) {
                    if (m < dd1MarkerOffset) {
                        ddNew.dosageValues[m][s] = dd1.dosageValues[m][s];
                    } else {
                        if (duplicatedSamplesAndIndices.containsKey(dd1.ids[s][0] + "\t" + dd1.ids[s][1])) {
                            ddNew.dosageValues[m][s] = dd2.dosageValues[dd2markersAndIndices.get(mkr)][duplicatedSamplesAndIndices.get(dd1.ids[s][0] + "\t" + dd1.ids[s][1])];
                        } else {
                            ddNew.dosageValues[m][s] = missingDosage;
                        }
                    }
                }
                int newInd = dd1.ids.length;
                for (int s = 0; s < dd2.ids.length; s++) {
                    if (duplicatedSampleIndices.contains(s)) {
                        continue;
                    }
                    if (!dd2markersAndIndices.containsKey(mkr)) {
                        ddNew.dosageValues[m][newInd] = missingDosage;
                    } else {
                        ddNew.dosageValues[m][newInd] = dd2.dosageValues[dd2IndM == -1 ? dd2markersAndIndices.get(mkr) : dd2IndM][s];
                    }
                    newInd++;
                }
                
                m++;
            }
        }
	    
	    return ddNew;
	}
	
	private DosageData() {
	    this.ids = null;
	    this.markerSet = null;
	    this.chrs = null;
	    this.positions = null;
	    this.alleles = null;
	    this.genotypeProbabilities = null;
	    this.dosageValues = null;
	}
	
	public DosageData(String dosageFile, String idFile, String mapFile, boolean verbose, Logger log) {
		this(dosageFile, idFile, mapFile, determineType(dosageFile), verbose, log);
	}
	
	public DosageData(String dosageFile, String idFile, String mapFile, int type, boolean verbose, Logger log) {
		this(dosageFile, idFile, mapFile, PARAMETERS[type], verbose, log);
	}
	
	public DosageData(String dosageFile, String idFile, String mapFile, int[] parameters, boolean verbose, Logger log) {
		BufferedReader reader;
		String[] line;
		String[] markerNames;
		Hashtable<String, String> invalids;
		
		markerSet = new SnpMarkerSet(mapFile);
		markerNames = markerSet.getMarkerNames();
		ids = HashVec.loadFileToStringMatrix(idFile, false, new int[] {0,1}, false);
		if (parameters[3] == 1) {
			dosageValues = new float[markerNames.length][ids.length];
		} else {
			genotypeProbabilities = new float[markerNames.length][ids.length][parameters[3]];
		}
		if (parameters[6] != -1) {
			alleles = new char[markerNames.length][2];
		}
		if (parameters[8] == CHR_INFO_IN_FILENAME) {
            Matcher m = Pattern.compile(CHR_REGEX).matcher(dosageFile);
            byte chr = -1;
            if (m.matches()) {
                chr = (byte) Integer.parseInt(m.group(1));
                String msg = "Warning - the format given expects chromosome number to be part of the file name.  This was determined to be chr{" + chr + "}.";
                if (log != null) {
                    log.report(msg);
                } else {
                    System.out.println(msg);
                }
                chrs = Array.byteArray(markerNames.length, chr);
            } else {
                String msg = "Error - the format given expects chromosome number to be part of the file name, but no chromosome number was found.  Chromosome information will not be included.";
                if (log != null) {
                    log.reportError(msg);
                } else {
                    System.err.println(msg);
                }
            }
		} else if (parameters[8] != -1) {
			chrs = new byte[markerNames.length];
		}
		if (parameters[9] != -1) {
			positions = new int[markerNames.length];
		}

		invalids = new Hashtable<String, String>();
		try {
			reader = Files.getAppropriateReader(dosageFile);//new BufferedReader(new FileReader(dosageFile));
			if (parameters[4] == 1) {
				line = reader.readLine().trim().split(parameters[10]==1?",":"[\\s]+");
				if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
					for (int i = 0; i < markerNames.length; i++) {
						if (!markerNames[i].equals(line[parameters[1]+parameters[3]*i]))	{
							String msg = "Error - mismatched name at marker "+(i+1)+" of "+dosageFile+"; expecting "+markerNames[i]+" given map file "+mapFile+", found "+line[parameters[1]+parameters[3]*i];
							if (log != null) {
							    log.reportError(msg);
							} else {
							    System.err.println(msg);
							}
							reader.close();
							return;
						}
					}
				} else if (parameters[2] == MARKER_DOMINANT_FORMAT) {
					if (parameters[3] != 2) {
					    String msg = "Warning - ignoring the header with IDs in file "+dosageFile+" because it does not contain 2 columns for each individual";
					    if (log != null) {
					        log.reportError(msg);
					    } else {
					        System.err.println(msg);
					    }
					} else {
						for (int i = 0; i < ids.length; i++) {
							if (!ids[i][0].equals(line[parameters[1]+parameters[3]*i+0]) || !ids[i][1].equals(line[parameters[1]+parameters[3]*i+1]))	{
							    String msg = "Error - mismatched IDs at individual "+(i+1)+" of "+dosageFile+"; expecting "+ids[i][0]+","+ids[i][1]+" given id file "+idFile+", found "+line[parameters[1]+parameters[3]*i+0]+","+line[parameters[1]+parameters[3]*i+1];
							    if (log != null) {
							        log.reportError(msg);
							    } else {
							        System.err.println(msg);
							    }
								reader.close();
								return;
							}
						}
					}
				}
			}

			if (parameters[2] == MARKER_DOMINANT_FORMAT) {
				for (int i = 0; i < markerNames.length; i++) {
					line = reader.readLine().trim().split(parameters[10]==1?",":"[\\s]+");
					if (!markerNames[i].equals(line[parameters[5]])) {
						String msg = "Error - mismatched name at marker "+(i+1)+" of "+dosageFile+"; expecting "+markerNames[i]+" given map file "+mapFile+", found "+line[parameters[5]];
						if (log != null) {
						    log.reportError(msg);
						} else {
						    System.err.println(msg);
						}
						reader.close();
						return;
					}

					if (parameters[6] != -1) {
						if (line[parameters[6]].length() > 1 || !Sequence.validAllele(line[parameters[6]])) {
							String msg = "Warning - invalid allele ('"+line[parameters[6]]+"') at marker "+markerNames[i];
	                        if (log != null) {
	                            log.reportError(msg);
	                        } else {
	                            System.err.println(msg);
	                        }
						}
						alleles[i][0] = line[parameters[6]].charAt(0);
					}
					if (parameters[7] != -1) {
						if (line[parameters[7]].length() > 1 || !Sequence.validAllele(line[parameters[7]])) {
							String msg = "Warning - invalid allele ('"+line[parameters[7]]+"') at marker "+markerNames[i];
	                        if (log != null) {
	                            log.reportError(msg);
	                        } else {
	                            System.err.println(msg);
	                        }
						}
						alleles[i][1] = line[parameters[7]].charAt(0);
					}
					if (parameters[8] >= 0) {
						try {
							chrs[i] = Byte.parseByte(line[parameters[8]]);
						} catch (NumberFormatException nfe) {
							chrs[i] = -1;
							if (!invalids.containsKey(line[parameters[8]])) {
								String msg = "Warning - invalid chromosome number ('"+line[parameters[8]]+"'), first seen at marker "+markerNames[i];
		                        if (log != null) {
		                            log.reportError(msg);
		                        } else {
		                            System.err.println(msg);
		                        }
								invalids.put(line[parameters[8]], ""); 
							}
						}
					}
					if (parameters[9] != -1) {
						try {
							positions[i] = Integer.parseInt(line[parameters[9]]);
						} catch (NumberFormatException nfe) {
							positions[i] = -1;
							if (!invalids.containsKey(line[parameters[9]])) {
								String msg = "Warning - invalid genome position ('"+line[parameters[9]]+"') for marker "+markerNames[i];
		                        if (log != null) {
		                            log.reportError(msg);
		                        } else {
		                            System.err.println(msg);
		                        }
								invalids.put(line[parameters[9]], ""); 
							}
						}
					}

					if (line.length - parameters[1] != ids.length * parameters[3]) {
					    String msg = "Error - mismatched number of elements in line "+(i+1+parameters[4])+" of "+dosageFile+"; expecting "+ids.length+"*"+parameters[3]+"+"+parameters[1]+"[="+(ids.length * parameters[3] + parameters[1])+"], found "+line.length;
                        if (log != null) {
                            log.reportError(msg);
                        } else {
                            System.err.println(msg);
                        }
						System.exit(1);
					}
					if (parameters[3] == 1) {
						for (int j = 0; j < ids.length; j++) {
							dosageValues[i][j] = Float.parseFloat(line[parameters[1]+j]);
						}
					} else {
						for (int j = 0; j < ids.length; j++) {
							genotypeProbabilities[i][j][0] = Float.parseFloat(line[parameters[1]+j*parameters[3]+0]);
							genotypeProbabilities[i][j][1] = Float.parseFloat(line[parameters[1]+j*parameters[3]+1]);
							if (parameters[3] == 3) {
								genotypeProbabilities[i][j][2] = Float.parseFloat(line[parameters[1]+j*parameters[3]+2]);
								if (Math.abs((1 - genotypeProbabilities[i][j][1] - genotypeProbabilities[i][j][0]) - genotypeProbabilities[i][j][2]) > 0.01) {
									String msg = "Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "+ids[j][0]+","+ids[j][1]+" at marker "+markerNames[i]+" which is line "+(i+1+parameters[4])+" of "+dosageFile+": "+line[parameters[1]+j*parameters[3]+0]+" "+line[parameters[1]+j*parameters[3]+1]+" "+line[parameters[1]+j*parameters[3]+2];
			                        if (log != null) {
			                            log.reportError(msg);
			                        } else {
			                            System.err.println(msg);
			                        }
								}
							}
						}
					}
				}
			} else if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
				for (int i = 0; i < ids.length; i++) {
					line = reader.readLine().trim().split(parameters[10]==1?",":"[\\s]+");
					if (line.length - parameters[1] != markerNames.length * parameters[3]) {
						String msg = "Error - mismatched number of elements in line "+(i+1+parameters[4])+" of "+dosageFile+"; expecting "+markerNames.length+"*"+parameters[3]+"+"+parameters[1]+"[="+(markerNames.length * parameters[3] + parameters[1])+"], found "+line.length;
                        if (log != null) {
                            log.reportError(msg);
                        } else {
                            System.err.println(msg);
                        }
						System.exit(1);
					}
					if (parameters[3] == 1) {
						for (int j = 0; j < markerNames.length; j++) {
							dosageValues[j][i] = Float.parseFloat(line[parameters[1]+j]);
						}
					} else {
						for (int j = 0; j < markerNames.length; j++) {
							genotypeProbabilities[j][i][0] = Float.parseFloat(line[parameters[1]+j*parameters[3]+0]);
							genotypeProbabilities[j][i][1] = Float.parseFloat(line[parameters[1]+j*parameters[3]+1]);
							if (parameters[3] == 3) {
								genotypeProbabilities[j][i][2] = Float.parseFloat(line[parameters[1]+j*parameters[3]+2]);
								if (Math.abs((1 - genotypeProbabilities[j][i][1] - genotypeProbabilities[j][i][0]) - genotypeProbabilities[j][i][2]) > 0.01) {
									String msg = "Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "+ids[i][0]+","+ids[i][1]+" at marker "+markerNames[j]+" which is line "+(i+1+parameters[4])+" of "+dosageFile+": "+line[parameters[1]+j*parameters[3]+0]+" "+line[parameters[1]+j*parameters[3]+1]+" "+line[parameters[1]+j*parameters[3]+2];
			                        if (log != null) {
			                            log.reportError(msg);
			                        } else {
			                            System.err.println(msg);
			                        }
								}
							}
						}
					}
				}
			}

			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dosageFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dosageFile + "\"");
			System.exit(2);
		}
	}
	
	public SnpMarkerSet getMarkerSet() {
		return markerSet;
	}

	public String[][] getIds() {
		return ids;
	}

	public float[][] getDosageValues() {
		return dosageValues;
	}

	public void computeDosageValues(Logger log) {
		if (genotypeProbabilities == null) {
			log.reportError("Error - cannot compute dosage values from genotype probabilities, if there are no genotype probabilities!");
			System.exit(1);
		}
		dosageValues = new float[genotypeProbabilities.length][genotypeProbabilities[0].length];
		for (int i = 0; i < dosageValues.length; i++) {
			for (int j = 0; j < ids.length; j++) {
				dosageValues[i][j] = genotypeProbabilities[i][j][0]*2+genotypeProbabilities[i][j][1]*1;
			}
		}
	}
	
	public void writeToFile(String filename, String mapOut, String extract, int format, Logger log) {
//		System.err.println("Error - type:" +format);
		writeToFile(filename, mapOut, extract, false, PARAMETERS[format], log);
	}	
	
	public void writeToFile(String filename, String mapOut, String extract, boolean allowIncompleteList, int[] parameters, Logger log) {
		PrintWriter writer;
		String[] line, markerNames;
//		char[][] alleles;
//		byte[] chrs;
//		int[] positions;
		String delimiter;
		String[] markersToKeep;
		HashSet<String> keeps;
		String root;
		SnpMarkerSet newMarkerSet; 

		
		if (extract == null) {
			keeps = null;
			markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut));
		} else {
			markersToKeep = HashVec.loadFileToStringArray(extract, false, new int[] {0}, false);
			keeps = HashVec.loadToHashSet(markersToKeep);
			root = ext.rootOf(filename, false);
			if (mapOut == null) {
				mapOut = root+".map";
			}
			newMarkerSet = markerSet.trim(markersToKeep, allowIncompleteList, false, log);
			if (!allowIncompleteList && newMarkerSet == null) {
				log.reportError("Error - failed to find all of the markers listed in '"+extract+"' in the file '"+filename+"'");
				return;
			}
			newMarkerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut));
			try {
				writer = new PrintWriter(new FileWriter(root+".ids.fam"));
				for (int i = 0; i < ids.length; i++) {
					writer.println(ids[i][0]+"\t"+ids[i][1]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + root+".ids.fam");
				e.printStackTrace();
			}
		}
		
		markerNames = markerSet.getMarkerNames();
		if (alleles==null) {
			alleles = markerSet.getAlleles();
		}
		if (chrs==null) {
			chrs = markerSet.getChrs();
			positions = markerSet.getPositions();
		}
		delimiter = DELIMITERS[parameters[10]];

		if (parameters[3] == 1 && dosageValues == null) {
			computeDosageValues(log);			
		}
		
		try {
			writer = new PrintWriter(new FileWriter(filename));
			if (parameters[4] == 1) {
				writer.print(Array.toStr(HEADS[parameters[11]], delimiter));
				
				if (parameters[2] == MARKER_DOMINANT_FORMAT) {
					if (parameters[3] == 2 && parameters[0] == FID_IID_TYPE) {
						for (int i = 0; i < ids.length; i++) {
							writer.print(delimiter+ids[i][0]+delimiter+ids[i][1]);
						}
					} else if (parameters[3] == 1 && parameters[0] == IID_TYPE) {
						for (int i = 0; i < ids.length; i++) {
							writer.print(delimiter+ids[i][1]);
						}						
					} else {
						log.reportError("Error - don't know how to list IDs when there "+(parameters[3]==1?"is one column":"are "+parameters[3]+" columns")+" for dosage inforation and the ID type is '"+parameters[2]+"'");
						System.exit(1);
					}
				} else if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
					for (int i = 0; i < markerNames.length; i++) {
						if (extract == null || keeps.contains(markerNames[i])) {
							for (int j = 0; j < parameters[3]; j++) {
								writer.print(delimiter+markerNames[i]);
							}
						}
					}
				}
				writer.println(); 
			}
			
			if (parameters[2] == MARKER_DOMINANT_FORMAT) {
				for (int i = 0; i < markerNames.length; i++) {
					if (extract == null || keeps.contains(markerNames[i])) {
						line = LEADS[parameters[11]]==null?new String[parameters[1]]:LEADS[parameters[11]];
						line[parameters[5]] = markerNames[i];
						if (parameters[6] >= 0) {
							if (alleles == null) {
								log.reportError("Error - file format requires alleles and none were supplied via the map file or the source dosage file");
								System.exit(1);
							}
							line[parameters[6]] = alleles[i][0]+"";
						}
						if (parameters[7] >= 0) {
							line[parameters[7]] = alleles[i][1]+"";
						}
						if (parameters[8] >= 0) {
							if (chrs == null) {
								log.reportError("Error - file format requires chrs/positions and none were supplied via the map file or the source dosage file");
								System.exit(1);
							}
							line[parameters[8]] = chrs[i]==-1?"-":chrs[i]+"";
						}
						if (parameters[9] >= 0) {
							line[parameters[9]] = positions[i]+"";
						}
						writer.print(Array.toStr(line, delimiter));
						
						if (parameters[3] == 1) {
							for (int j = 0; j < ids.length; j++) {
								writer.print(delimiter+ext.formDeci(dosageValues[i][j], parameters[13], parameters[12] == parameters[13]));
							}
						} else {
							if (genotypeProbabilities == null) {
								log.reportError("Error - cannot write genotype probabilities if they were never read in!");
								System.exit(1);
							}
							for (int j = 0; j < ids.length; j++) {
								writer.print(delimiter+ext.formDeci(genotypeProbabilities[i][j][0], parameters[13], parameters[12] == parameters[13]));
								writer.print(delimiter+ext.formDeci(genotypeProbabilities[i][j][1], parameters[13], parameters[12] == parameters[13]));
								if (parameters[3] == 3) {
									if (genotypeProbabilities[i][j].length > 2) {
										writer.print(delimiter+ext.formDeci(genotypeProbabilities[i][j][2], parameters[13], parameters[12] == parameters[13]));
									} else {
										writer.print(delimiter+ext.formDeci(1 - genotypeProbabilities[i][j][1] - genotypeProbabilities[i][j][0], parameters[13], parameters[12] == parameters[13]));
									}
								}
							}
						}
						writer.println();
					}
				}
			} else if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
				for (int i = 0; i < ids.length; i++) {
					line = LEADS[parameters[11]]==null?new String[parameters[3]]:LEADS[parameters[11]];
					if (parameters[0] == MACH_ID_TYPE) {
						line[parameters[5]] = ids[i][0]+"->"+ids[i][1];
					} else if (parameters[0] == IID_TYPE) {
						line[parameters[5]] = ids[i][1];
					} else if (parameters[0] == FID_IID_TYPE) {
						line[parameters[5]] = ids[i][0]+delimiter+ids[i][1];
					} else {
						log.reportError("Error - ID type has not been defined");
						System.exit(1);
					}
					writer.print(Array.toStr(line, delimiter));
					                                  
					if (parameters[3] == 1) {
						for (int j = 0; j < markerNames.length; j++) {
							if (extract == null || keeps.contains(markerNames[j])) {
								writer.print(delimiter+ext.formDeci(dosageValues[j][i], parameters[13], parameters[12] == parameters[13]));
							}
						}
					} else {
						for (int j = 0; j < markerNames.length; j++) {
							if (extract == null || keeps.contains(markerNames[j])) {
								writer.print(delimiter+ext.formDeci(genotypeProbabilities[j][i][0], parameters[13], parameters[12] == parameters[13]));
								writer.print(delimiter+ext.formDeci(genotypeProbabilities[j][i][1], parameters[13], parameters[12] == parameters[13]));
								if (parameters[3] == 3) {
									if (genotypeProbabilities[j][i].length > 2) {
										writer.print(delimiter+ext.formDeci(genotypeProbabilities[j][i][2], parameters[13], parameters[12] == parameters[13]));
									} else {
										writer.print(delimiter+ext.formDeci(1 - genotypeProbabilities[j][i][1] - genotypeProbabilities[j][i][0], parameters[13], parameters[12] == parameters[13]));
									}
								}
							}
						}
					}
					writer.println();
				}
			}
			
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public void analyze(String phenoFile, String phenoMissingValue, String snpList, boolean verbose, Logger log) {
		PrintWriter writer, w2;
		String[] line;
		Hashtable<String, String> hash;
		HashSet<String> snps;
		int count;
		String[] traits, markerNames;
		boolean[] use, analyze;
		byte[] chrs;
		int[] positions;
		double[] deps;
		double[][] indeps;
		boolean logistic;
		RegressionModel model;
		char[][] alleles;
		double[] betas, stderrs, pvals, stats;

		traits = Files.getHeaderOfFile(phenoFile, "\t", log);
		hash = HashVec.loadFileToHashString(phenoFile, new int[] {0,1}, Array.subArray(Array.intArray(traits.length), 2, traits.length), false, "\t", true, false, false);
		traits = Array.subArray(traits, 2);

		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		alleles = markerSet.getAlleles();
		analyze = Array.booleanArray(markerNames.length, true);
		if (snpList != null) {
			snps = HashVec.loadFileToHashSet(snpList, false);
			for (int i = 0; i < markerNames.length; i++) {
				if (!snps.contains(markerNames[i])) {
					analyze[i] = false;
				}
			}
		}
		
		use = Array.booleanArray(ids.length, true);
		for (int i = 0; i < ids.length; i++) {
			if (hash.containsKey(ids[i][0]+"\t"+ids[i][1])) {
				line = hash.get(ids[i][0]+"\t"+ids[i][1]).split("[\\s]+");
				for (int j = 0; j < line.length; j++) {
					if (!ext.isValidDouble(line[j]) && !line[j].equals(phenoMissingValue)) {
						use[i] = false;
					}
				}
			} else {
				use[i] = false;
			}
		}
		deps = new double[Array.booleanArraySum(use)];
		indeps = new double[deps.length][traits.length];
		log.report("There are "+deps.length+" rows with complete data", true, verbose);
		count = 0;
		for (int i = 0; i < ids.length; i++) {
			if (use[i]) {
				line = hash.get(ids[i][0]+"\t"+ids[i][1]).split("[\\s]+");
				deps[count] = Double.parseDouble(line[0]);
				for (int j = 1; j < traits.length; j++) {
					indeps[count][j] = Double.parseDouble(line[j]);
				}
				count++;
			}
		}
		logistic = RegressionModel.isBinaryTrait(Array.toStr(deps).split("[\\s]+"), log);
		log.report("Running a "+(logistic?"logistic":"linear")+" model for trait '"+traits[0]+"'", true, verbose);
		try {
//			writer = new PrintWriter(new FileWriter(ext.rootOf(phenoFile, false)+(snpList==null?"":"_"+ext.rootOf(snpList, false))+".results."+(logistic?"logistic":"linear")));
//			w2 = new PrintWriter(new FileWriter(ext.rootOf(phenoFile, false)+(snpList==null?"":"_"+ext.rootOf(snpList, false))+".se.metal"));
			writer = new PrintWriter(new FileWriter(ext.rootOf(phenoFile, false)+".results."+(logistic?"logistic":"linear")));
			w2 = new PrintWriter(new FileWriter(ext.rootOf(phenoFile, false)+".se.metal"));
			line = Array.clone(logistic?Plink.LOGISTIC_SE_HEADER:Plink.LINEAR_SE_HEADER);
			line[1] = line[1]+"      ";
			line[2] = line[1]+"      ";
			writer.println(Array.toStr(line));
//			w2.println("MARKER\tREF\tOTHER\tN\tDIR\tPVALUE\tbeta\tSE");
			w2.println("MarkerName\tAllele1\tAllele2\tWeight\tDirection\tP-value\tEffect\tStdErr");
//			public static final String[] LOGISTIC_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "SE", "L95", "U95", "STAT", "P"};
			for (int i = 0; i < markerNames.length; i++) {
				if (analyze[i]) {
					count = 0;
					for (int j = 0; j < ids.length; j++) {
						if (use[j]) {
							indeps[count][0] = (double)dosageValues[i][j];
							count++;
						}
					}
					model = logistic?new LogisticRegression(deps, indeps, false, false):new LeastSquares(deps, indeps, false, false);
					betas = model.getBetas();
					stderrs = model.getSEofBs();
					pvals = model.getSigs();
					stats = model.getStats();
					int sigfig = 4;
//					System.err.println(betas.length+"\t"+traits.length);
					if (betas.length != traits.length+1) {
						writer.println(chrs[i]+"\t"+markerNames[i]+"\t"+positions[i]+"\t"+alleles[i][0]+"\tADD\t.\t.\t.\t.\t.\t.\t.");
						w2.println(markerNames[i]+"\t"+alleles[i][0]+"\t"+alleles[i][1]+"\t"+deps.length+"\t.\t.\t.\t.");
					} else {
						writer.println(chrs[i]+"\t"+markerNames[i]+"\t"+positions[i]+"\t"+alleles[i][0]+"\tADD\t"+deps.length+"\t"+ext.formDeci(betas[1], sigfig, true)+"\t"+ext.formDeci(stderrs[1], sigfig, true)+"\t.\t.\t"+(stats[1]+"    ").substring(0,6).trim()+"\t"+ext.prettyP(pvals[1], sigfig, 4, 3, true));
						w2.println(markerNames[i]+"\t"+alleles[i][0]+"\t"+alleles[i][1]+"\t"+deps.length+"\t"+(betas[1]==0?0:(betas[1]>0?"+":"-"))+"\t"+ext.prettyP(pvals[1], sigfig, 4, 3, true)+"\t"+ext.formDeci(betas[1], 6, true)+"\t"+ext.formDeci(stderrs[1], 6, true));
						for (int j = 1; j < traits.length; j++) {
							writer.println(chrs[i]+"\t"+markerNames[i]+"\t"+positions[i]+"\t"+alleles[i][0]+"\t"+traits[j]+"\t"+deps.length+"\t"+ext.formDeci(betas[1+j], sigfig, true)+"\t"+ext.formDeci(stderrs[1+j], sigfig, true)+"\t.\t.\t"+(stats[1+j]+"     ").substring(0,6).trim()+"\t"+ext.prettyP(pvals[1+j], sigfig, 4, 3, true));
						}
					}
					writer.flush();
				}
			}
			writer.close();
			w2.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(phenoFile, false)+(snpList==null?"":"_"+ext.rootOf(snpList, false))+".results."+(logistic?"logistic":"linear"));
			e.printStackTrace();
		}
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static DosageData load(String filename, boolean jar) {
		return (DosageData)Files.readSerial(filename, jar, true);
	}

	public static int determineType(String dosageFileName) {
		String dosageFile = dosageFileName;
		if (dosageFile.endsWith(".gz")) {
			dosageFile = dosageFile.substring(0, dosageFile.length()-3);
		} else if (dosageFile.endsWith(".zip")) {
			dosageFile = dosageFile.substring(0, dosageFile.length()-4);
		}
		if (dosageFile.endsWith(".mldose")) {
			return MACH_MLDOSE_FORMAT;
		} else if (dosageFile.endsWith(".gen")) {
			return GEN_FORMAT;
		} else if (dosageFile.endsWith(".fhsR")) {
			return GWAF_FORMAT;
		} else if (dosageFile.endsWith(".dosage")) {
			return PLINK_FORMAT;
		} else if (dosageFile.endsWith(".mlprob")) {
			return MACH_MLPROB_FORMAT;
		} else if (dosageFile.endsWith(".dose")) {
			return MINIMAC_DOSE_FORMAT;
		} else if (dosageFile.endsWith(".impute2") || dosageFile.endsWith(".imputed")) {
			return IMPUTE2_DOSE_FORMAT;
		} else if (dosageFile.endsWith(".db.xln")) {
			return DATABASE_DOSE_FORMAT;
		} else {
			System.err.println("Error - format could not be deduced solely by the filename extension");
			return -1;
		}
	}
	
	public static void convert(String dosageFile, String idFile, String mapFile, int fromFormat, String outfile, String mapOut, String extract, int toFormat, boolean awk, boolean verbose, Logger log) {
		convert(dosageFile, idFile, mapFile, PARAMETERS[fromFormat], outfile, mapOut, extract, PARAMETERS[toFormat], awk, verbose, log);
	}
	
	public static void convert(String dosageFile, String idFile, String mapFile, int[] fromParameters, String outfile, String mapOut, String extract, int[] toParameters, boolean awk, boolean verbose, Logger log) {
		SnpMarkerSet markerSet;
		String[][] ids;
		float dosageValue;
		float[] genotypeProbabilities;
		byte[] chrs;
		int[] positions;
		char[][] alleles;
		byte chr;
		int position;
		char[] allelePair;
		
		BufferedReader reader;
		String[] line, markerNames;
		Hashtable<String, String> invalids;

		PrintWriter writer;
		String delimiter;
		String[] lead, markersToKeep;
		HashSet<String> keeps;
		String root;
		IntVector cols;
		int offset, numColsEach;
		boolean dos;
		String ranges;
		int start, stop;
		
		
		if (fromParameters[2] != toParameters[2]) {
			log.reportError("Conversion will have to take place in memory in order to transpose between marker dominant/individual dominant; use new DosageDate().writeToFile() instead");
			return;
		}
		
		markerSet = new SnpMarkerSet(mapFile, false, new Logger());
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		alleles = markerSet.getAlleles();
		ids = HashVec.loadFileToStringMatrix(idFile, false, new int[] {0,1}, false);
		
		if (extract == null) {
			keeps = null;
		} else {
			markersToKeep = HashVec.loadFileToStringArray(extract, false, new int[] {0}, false);
			keeps = HashVec.loadToHashSet(markersToKeep);
			root = ext.rootOf(outfile, false);
			markerSet = markerSet.trim(markersToKeep, true, false, log);		// allows missing markers, but will list how many
			if (mapOut == null) {
				mapOut = root+".pinfo";
				markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut));
				mapOut = root+".map";
			}
			markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut));
			try {
				writer = new PrintWriter(new FileWriter(root+".ids.fam"));
				for (int i = 0; i < ids.length; i++) {
					writer.println(ids[i][0]+"\t"+ids[i][1]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + root+".ids.fam");
				e.printStackTrace();
			}
		}
		
//		/** 0,       1,                                                  2,                3,                                      4,          5,                6,        7,        8,         9,         10,        11,                  12,                   13                   */
//		/** id_type, column index where dosage/probability values begin, dominance format, number of columns summarizing the data, header row, marker/IID index, A1 index, A2 index, chr index, pos index, delimiter, index for head/lead, min number of digits, max number of digits */
//
//		/** 									   0                      1  2                           3  4  5   6   7   8   9 10 11 12 13 */
//		public static final int[][] PARAMETERS = {{MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1, -1, -1, -1, 0, 0, 3, 3}, // .mldose (MACH)
//												  {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT,     3, 0, 1,  3,  4,  0,  2, 2, 1, 0, 3}, // .gen
//												  {IID_TYPE,              1, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0, -1, -1, -1, -1, 1, 2, 0, 3}, // .fhsR (GWAF) 
//												  {FID_IID_TYPE,          3, MARKER_DOMINANT_FORMAT,     2, 1, 0,  1,  2, -1, -1, 0, 3, 2, 2}, // .dosage (PLINK)
//												  {MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 2, 0, 0, -1, -1, -1, -1, 0, 4, 3, 3}, // .mlprob (MACH)
//												  {MACH_ID_TYPE,          2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1, -1, -1, -1, 0, 5, 3, 3}, // .dose (MINIMAC)
//		};
		
		
		if (awk && !Array.equals(fromParameters, toParameters)) {
			log.reportError("Error - the awk option of convert is currently only available for identical file types");
		} else if (awk) {
			cols = new IntVector();
			offset = fromParameters[1];
			for (int i = 2; i <= offset; i++) {
				cols.add(i);
			}
			numColsEach = fromParameters[3];
			for (int i = 0; i < markerNames.length; i++) {
				if (keeps.contains(markerNames[i])) {
					for (int j = 0; j < numColsEach; j++) {
						cols.add(offset+i*numColsEach+j+1);
					}
				}
			}
			ranges = ext.listRanges(cols.toArray());
			System.out.println("Extracting columns: 1,"+ranges);
			
			dos = System.getProperty("os.name").startsWith("Windows");
			try {
				writer = new PrintWriter(new FileWriter("awk_command.bat"));
				writer.print("awk "+(dos?"\"":"'")+"BEGIN { OFS = "+(dos?"\\":"")+"\"\\t"+(dos?"\\":"")+"\" } { printf $1 ;");
				line = ranges.split(",");
				for (int i = 0; i < line.length; i++) {
					if (line[i].contains("-")) {
						start = Integer.parseInt(line[i].substring(0, line[i].indexOf("-")));
						stop = Integer.parseInt(line[i].substring(line[i].indexOf("-")+1));
					} else {
						start = stop = Integer.parseInt(line[i]);
					}
					writer.print(" for (i="+start+"; i<="+stop+"; i++) printf "+(dos?"\\":"")+"\"\\t"+(dos?"\\":"")+"\"$i ;");
				}
				writer.println(" printf "+(dos?"\\":"")+"\"\\n"+(dos?"\\":"")+"\" }"+(dos?"\"":"'")+" "+dosageFile +" > "+outfile);
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + "awk_command.bat");
				e.printStackTrace();
			}

			
//			System.out.println("Running the following command to parse the new dosage file:\n"+command);
			Files.chmod("awk_command.bat", false);
			CmdLine.run((dos?"":"./")+"awk_command.bat", "./", System.err);
		} else {
			allelePair = new char[2];
			if (fromParameters[3] > 1) {
				genotypeProbabilities = new float[fromParameters[3]];
			} else {
				genotypeProbabilities = null;
			}
			dosageValue = -999;
			
			delimiter = DELIMITERS[toParameters[10]];
			invalids = new Hashtable<String, String>();
			try {
				reader = Files.getAppropriateReader(dosageFile); // new BufferedReader(new FileReader(dosageFile));
				writer = new PrintWriter(new FileWriter(outfile));
	
				if (fromParameters[4] == 1) {
					line = reader.readLine().trim().split(fromParameters[10]==1?",":"[\\s]+");
					if (fromParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
						for (int i = 0; i < markerNames.length; i++) {
							if (!markerNames[i].equals(line[fromParameters[1]+fromParameters[3]*i]))	{
								log.reportError("Error - mismatched name at marker "+(i+1)+" of "+dosageFile+"; expecting "+markerNames[i]+" given map file "+mapFile+", found "+line[fromParameters[1]+fromParameters[3]*i]);
								reader.close();
								return;
							}
						}
					} else if (fromParameters[2] == MARKER_DOMINANT_FORMAT) {
						if (fromParameters[3] != 2) {
							log.reportError("Warning - ignoring the header with IDs in file "+dosageFile+" because it does not contain 2 columns for each individual");
						} else {
							for (int i = 0; i < ids.length; i++) {
								if (!ids[i][0].equals(line[fromParameters[1]+fromParameters[3]*i+0]) || !ids[i][1].equals(line[fromParameters[1]+fromParameters[3]*i+1]))	{
									log.reportError("Error - mismatched IDs at individual "+(i+1)+" of "+dosageFile+"; expecting "+ids[i][0]+","+ids[i][1]+" given id file "+idFile+", found "+line[fromParameters[1]+fromParameters[3]*i+0]+","+line[fromParameters[1]+fromParameters[3]*i+1]);
									reader.close();
									return;
								}
							}
						}
					}
				}
	
				if (toParameters[4] == 1) { // if there is a header row
					writer.print(Array.toStr(HEADS[toParameters[11]], delimiter));
					
					if (toParameters[2] == MARKER_DOMINANT_FORMAT) {
						if (toParameters[3] == 2 && toParameters[0] == FID_IID_TYPE) {
							for (int i = 0; i < ids.length; i++) {
								writer.print(delimiter+ids[i][0]+delimiter+ids[i][1]);
							}
						} else if (toParameters[3] == 1 && toParameters[0] == IID_TYPE) {
							for (int i = 0; i < ids.length; i++) {
								writer.print(delimiter+ids[i][1]);
							}						
						} else {
							log.reportError("Error - don't know how to list IDs when there "+(toParameters[3]==1?"is one column":"are "+toParameters[3]+" columns")+" for dosage inforation and the ID type is '"+toParameters[2]+"'");
							System.exit(1);
						}
					} else if (toParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
						for (int i = 0; i < markerNames.length; i++) {
							if (extract == null || keeps.contains(markerNames[i])) {
								for (int j = 0; j < toParameters[3]; j++) {
									writer.print(delimiter+markerNames[i]);
								}
							}
						}
					}
					writer.println(); 
				}
	
				if (fromParameters[2] == MARKER_DOMINANT_FORMAT) {
					for (int i = 0; i < markerNames.length; i++) {
						line = reader.readLine().trim().split(fromParameters[10]==1?",":"[\\s]+");
						if (!markerNames[i].equals(line[fromParameters[5]])) {
							log.reportError("Error - mismatched name at marker "+(i+1)+" of "+dosageFile+"; expecting "+markerNames[i]+" given map file "+mapFile+", found "+line[fromParameters[5]]);
							reader.close();
							return;
						}
	
						// gathering alleles if available
						if (fromParameters[6] != -1) {
							if (line[fromParameters[6]].length() > 1 || !Sequence.validAllele(line[fromParameters[6]])) {
								log.reportError("Warning - invalid allele ('"+line[fromParameters[6]]+"') at marker "+markerNames[i]);
							}
							allelePair[0] = line[fromParameters[6]].charAt(0);
						} else if (alleles != null) {
							allelePair[0] = alleles[i][0];
						} else {
							allelePair[0] = '~';
						}
	
						// gathering alleles if available
						if (fromParameters[7] != -1) {
							if (line[fromParameters[7]].length() > 1 || !Sequence.validAllele(line[fromParameters[7]])) {
								log.reportError("Warning - invalid allele ('"+line[fromParameters[7]]+"') at marker "+markerNames[i]);
							}
							allelePair[1] = line[fromParameters[7]].charAt(0);
						} else if (alleles != null) {
							allelePair[1] = alleles[i][1];
						} else {
							allelePair[1] = '~';
						}
	
						// gathering chromosome if available
						if (fromParameters[8] != -1) {
							try {
								chr = Byte.parseByte(line[fromParameters[8]]);
							} catch (NumberFormatException nfe) {
								chr = -1;
								if (!invalids.containsKey(line[fromParameters[8]])) {
									log.reportError("Warning - invalid chromosome number ('"+line[fromParameters[8]]+"'), first seen at marker "+markerNames[i]);
									invalids.put(line[fromParameters[8]], ""); 
								}
							}
						} else if (chrs != null) {
							chr = chrs[i];
						} else {
							chr = -9;
						}
	
						// gathering position if available
						if (fromParameters[9] != -1) {
							try {
								position = Integer.parseInt(line[fromParameters[9]]);
							} catch (NumberFormatException nfe) {
								position = -1;
								if (!invalids.containsKey(line[fromParameters[9]])) {
									log.reportError("Warning - invalid genome position ('"+line[fromParameters[9]]+"') for marker "+markerNames[i]);
									invalids.put(line[fromParameters[9]], ""); 
								}
							}
						} else if (positions != null) {
							position = positions[i];
						} else {
							position = 0;
						}
	
						if (line.length - fromParameters[1] != ids.length * fromParameters[3]) {
							log.reportError("Error - mismatched number of elements in line "+(i+1+fromParameters[4])+" of "+dosageFile+"; expecting "+ids.length+"*"+fromParameters[3]+"+"+fromParameters[1]+", found "+line.length);
							System.exit(1);
						}
						
						if (extract == null || keeps.contains(markerNames[i])) {
							lead = LEADS[toParameters[11]]==null?new String[toParameters[1]]:LEADS[toParameters[11]];
							lead[toParameters[5]] = markerNames[i];

							// adding alleles if required
							if (toParameters[6] >= 0) {
								if (allelePair[0] == '~' || allelePair[1] == '~') {
									log.reportError("Error - file format requires alleles and none were supplied via the map file or the source dosage file");
									System.exit(1);
								}
								lead[toParameters[6]] = allelePair[0]+"";
							}
							// adding alleles if required
							if (toParameters[7] >= 0) {
								lead[toParameters[7]] = allelePair[1]+"";
							}
							// adding chromosome if required
							if (toParameters[8] >= 0) {
								if (chr == -9) {
									log.reportError("Error - file format requires chrs/positions and none were supplied via the map file or the source dosage file");
									System.exit(1);
								}
								lead[toParameters[8]] = chr==-1?"-":chr+"";
							}
							// adding position if required
							if (toParameters[9] >= 0) {
								lead[toParameters[9]] = position+"";
							}
							writer.print(Array.toStr(lead, delimiter));
	
						
							for (int j = 0; j < ids.length; j++) {
								if (fromParameters[3] == 1) {
									dosageValue = Float.parseFloat(line[fromParameters[1]+j]);
								} else {
									genotypeProbabilities[0] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+0]);
									genotypeProbabilities[1] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+1]);
									if (fromParameters[3] == 3) {
										genotypeProbabilities[2] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+2]);
										if (Math.abs((1 - genotypeProbabilities[1] - genotypeProbabilities[0]) - genotypeProbabilities[2]) > 0.01) {
											log.reportError("Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "+ids[j][0]+","+ids[j][1]+" at marker "+markerNames[i]+" which is line "+(i+1+fromParameters[4])+" of "+dosageFile+": "+line[fromParameters[1]+j*fromParameters[3]+0]+" "+line[fromParameters[1]+j*fromParameters[3]+1]+" "+line[fromParameters[1]+j*fromParameters[3]+2]);
										}
									}
								}
								
								if (toParameters[3] == 1) {
									if (fromParameters[3] > 1) {
										dosageValue = genotypeProbabilities[0]*2+genotypeProbabilities[1]*1;
									}
									writer.print(delimiter+ext.formDeci(dosageValue, toParameters[13], toParameters[12] == toParameters[13]));
								} else {
									writer.print(delimiter+ext.formDeci(genotypeProbabilities[0], toParameters[13], toParameters[12] == toParameters[13]));
									writer.print(delimiter+ext.formDeci(genotypeProbabilities[1], toParameters[13], toParameters[12] == toParameters[13]));
									if (toParameters[3] == 3) {
										if (fromParameters[3] == 3) {
											writer.print(delimiter+ext.formDeci(genotypeProbabilities[2], toParameters[13], toParameters[12] == toParameters[13]));
										} else {
											writer.print(delimiter+ext.formDeci(1 - genotypeProbabilities[1] - genotypeProbabilities[0], toParameters[13], toParameters[12] == toParameters[13]));
										}
									}
								}
							}
							writer.println();
						}
					}
				} else if (fromParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
					for (int i = 0; i < ids.length; i++) {
						line = reader.readLine().trim().split(fromParameters[10]==1?",":"[\\s]+");
						if (line.length - fromParameters[1] != markerNames.length * fromParameters[3]) {
							log.reportError("Error - mismatched number of elements in line "+(i+1+fromParameters[4])+" of "+dosageFile+"; expecting "+markerNames.length+"*"+fromParameters[3]+"+"+fromParameters[1]+", found "+line.length);
							System.exit(1);
						}
						
						lead = LEADS[toParameters[11]]==null?new String[toParameters[3]]:LEADS[toParameters[11]];
						if (toParameters[0] == MACH_ID_TYPE) {
							lead[toParameters[5]] = ids[i][0]+"->"+ids[i][1];
						} else if (toParameters[0] == IID_TYPE) {
							lead[toParameters[5]] = ids[i][1];
						} else if (toParameters[0] == FID_IID_TYPE) {
							lead[toParameters[5]] = ids[i][0]+delimiter+ids[i][1];
						} else {
							log.reportError("Error - ID type has not been defined");
							System.exit(1);
						}
						writer.print(Array.toStr(lead, delimiter));
						
						for (int j = 0; j < markerNames.length; j++) {
							if (extract == null || keeps.contains(markerNames[j])) {
								if (fromParameters[3] == 1) {
									dosageValue = Float.parseFloat(line[fromParameters[1]+j]);
								} else {
									genotypeProbabilities[0] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+0]);
									genotypeProbabilities[1] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+1]);
									if (fromParameters[3] == 3) {
										genotypeProbabilities[2] = Float.parseFloat(line[fromParameters[1]+j*fromParameters[3]+2]);
										if (Math.abs((1 - genotypeProbabilities[1] - genotypeProbabilities[0]) - genotypeProbabilities[2]) > 0.01) {
											log.reportError("Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "+ids[i][0]+","+ids[i][1]+" at marker "+markerNames[j]+" which is line "+(i+1+fromParameters[4])+" of "+dosageFile+": "+line[fromParameters[1]+j*fromParameters[3]+0]+" "+line[fromParameters[1]+j*fromParameters[3]+1]+" "+line[fromParameters[1]+j*fromParameters[3]+2]);
										}
									}
								}
	
								if (toParameters[3] == 1) {
									if (fromParameters[3] > 1) {
										dosageValue = genotypeProbabilities[0]*2+genotypeProbabilities[1]*1;
									}
									writer.print(delimiter+ext.formDeci(dosageValue, toParameters[13], toParameters[12] == toParameters[13]));
								} else {
									writer.print(delimiter+ext.formDeci(genotypeProbabilities[0], toParameters[13], toParameters[12] == toParameters[13]));
									writer.print(delimiter+ext.formDeci(genotypeProbabilities[1], toParameters[13], toParameters[12] == toParameters[13]));
									if (toParameters[3] == 3) {
										if (genotypeProbabilities.length > 2) {
											writer.print(delimiter+ext.formDeci(genotypeProbabilities[2], toParameters[13], toParameters[12] == toParameters[13]));
										} else {
											writer.print(delimiter+ext.formDeci(1 - genotypeProbabilities[1] - genotypeProbabilities[0], toParameters[13], toParameters[12] == toParameters[13]));
										}
									}
								}
							}
						}
						writer.println();
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + dosageFile + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + dosageFile + "\"");
				System.exit(2);
			}
		}		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "chr21.mldose";
		String idFile = "chr21.fam";
		String mapFile = "chr21.mlinfo";
		String outfile = "output.dosage";
		String mapOut = null;
		int from = -1;
		int to = -1;
		String logfile = null;
		Logger log;
		String extract = null;
		boolean awk = false;
		long date;

		String usage = "\n" +
		"filesys.DosageData requires 0-1 arguments\n" +
		"   (1) name of dosage file to convert (i.e. dosageFile=" + filename + " (default))\n" + 
		"   (2) name of file with ids listed (i.e. idFile=" + idFile + " (default))\n" + 
		"   (3) name of associated map file (i.e. mapFile=" + mapFile + " (default))\n" + 
		"   (4) name of new dosage file (i.e. out=" + outfile + " (default))\n" + 
		"   (5) name of new map file (i.e. mapOut=" + mapOut + " (default))\n" + 
		"   (6) filename for log (i.e. log=" + logfile + " (default) (optional))\n" + 
		"   (7) name of file listing variants to extract (i.e. extract=" + extract + " (default) (optional))\n" + 
		"   (8) file type of original file (i.e. from=" + from + " (default))\n" + 
		"   (9) file type of file to be generated (i.e. to=" + to + " (default))\n" + 
		"  (10) use sed/cut instead (i.e. -awk (not the default))\n" + 
		" \n" + 
		"       Type Extension Description\n" + 
		"        -1   [any]    auto-detect from extension\n" + 
		"        0   .mldose   MACH style .mldose file\n" + 
		"        1   .gen      probabilites for all three genotypes, plus map information\n" + 
		"        2   .fhsR     input file for GWAF\n" + 
		"        3   .dosage   input file for PLINK\n" + 
		"        4   .mlprob   MACH style .mlprob file\n" + 
		"        5   .dose     output from MINIMAC\n" + 
		"        6   .impute2  output from IMPUTE2\n" + 
		"        7   .db.xln   database format\n" + 
		"        8   .dose     BEAGLE format\n" + 
		"        [Files may also be zipped or gzipped]" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dosageFile=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("idFile=")) {
				idFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("mapFile=")) {
				mapFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("mapOut=")) {
				mapOut = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("extract=")) {
				extract = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("from=")) {
				from = ext.parseIntArg(args[i]);				
				numArgs--;
			} else if (args[i].startsWith("to=")) {
				to = ext.parseIntArg(args[i]);				
				numArgs--;
			} else if (args[i].startsWith("-awk")) {
				awk = true;				
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
//		filename = "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\00src\\leslie_lange.CFS.IBC.CEU.chr9.gen";
//		idFile = "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\leslie_lange.CFS.IBC.CEU.chr9.pfam";
//		mapFile = "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\leslie_lange.CFS.IBC.CEU.chr9.pmap";
//		outfile = "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\00src\\leslie_lange.CFS.IBC.CEU.chr9.fhsR";

//		filename = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.gen";
//		idFile = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.pfam";
//		mapFile = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.mlinfo";
//		outfile = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\candi\\abo.gen";
//		extract = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\abo_snplist.txt";
		
//		filename = "chr4.dose";
//		mapFile = "chr4.minfo";
//		idFile = "plink.fam";
//		outfile = "first_half.dose";
//		extract = "first_half.txt";
//		awk = true;
//
//		filename = "Fung/hla.dose";
//		mapFile = "Fung/hla.minfo";
//		idFile = "Fung/hla.ids.fam";
//		outfile = "hla_test.dose";
//		extract = "test_list.txt";
//		awk = true;
		
//		filename = "D:/LITE/1000Gimputation/blacks/chr9.all_parsed.xln";
//		idFile = "D:/LITE/1000Gimputation/blacks/B_ARIC_gwas_frz3v2_clean.b37.fam";
//		mapFile = "D:/LITE/1000Gimputation/blacks/ABO_blacks_fullRegion.bim";
//		outfile = "D:/LITE/1000Gimputation/blacks/ABO_blacks_fullRegion.db.xln";
//		from=6;
//		to=7;

		from = from==-1?determineType(filename):from;
		to = to==-1?determineType(outfile):to;
		try {
			log = new Logger(logfile);
			if (PARAMETERS[from][3] < PARAMETERS[to][3]) {
				log.reportError("Error - cannot convert dosage values to genotype probabilities");
			} else if (PARAMETERS[from][2] != PARAMETERS[to][2]) {
				log.reportError("Conversion will have to take place in memory in order to transpose between marker dominant/individual dominant");
				new DosageData(filename, idFile, mapFile, from, true, log).writeToFile(outfile, mapOut, extract, to, log);
			} else {
				date = new Date().getTime();
				convert(filename, idFile, mapFile, from, outfile, mapOut, extract, to, awk, true, log);
				System.out.println("Time elapsed: "+ext.getTimeElapsed(date));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
