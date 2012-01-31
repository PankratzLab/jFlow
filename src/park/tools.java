package park;

import java.io.*;
import java.util.*;

import common.*;

public class tools {
	public static final String NINFO_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\00ninfos\\";
	public static final String CRF_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\";
	public static final String PARKIN_DIR = CRF_DIR+"parkin\\";
	public static final String MASTER_PLATELIST_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\2_platelists\\";
	public static final String[] NINFO1_HEADER = {"FamNo", "IndNo", "AgeOfOnset", "IllnessStatusCode", "DNA", "Sex", "IllnessCode", "RaceCode"};
	public static final String[] NINFO2_HEADER = {"FamNo", "IndNo", "Sex", "Deceased", "FatherInd", "MotherInd", "IllnessCode", "IllnessStatusCode"};
	public static final String[] PD_NAMES = {"PD", "PD+", "LB", "DLB"};
	public static final String[] PD_DX = {"CONF_PD", "VPD", "NVPD", "RPD"};
	public static final String[] AFFECTED = {"VPD", "NVPD", "RPD", "FRPT", "CONF", "POSS", "PROB"};
	public static final String[] UNAFFECTED = {"NOEV", "NRPD", "NXAM"};
	public static final String ALT_DIR = "/home/npankrat/park/00masters/crf_db/";
	public static final String DB_FILE = "crf_db.dat";

	public static String getUniqueID(String famid, String indid) {
		if (indid.startsWith("ND")) {
			try {
				Integer.parseInt(indid.substring(2));
				return indid;
			} catch (NumberFormatException e) {}
		}

		if (famid.length()!=5||indid.length()>3) {
			System.err.println("Warning - '"+famid+"' and '"+indid+"', wrong lengths for a PD identifier");
		}
		return famid+ext.formNum(Integer.parseInt(indid), 3);
	}

	public static String[] getFamID(String famIndId) {
		int index = famIndId.indexOf("-");

		if (index==-1&&famIndId.length()==8&&(famIndId.startsWith("7")||famIndId.startsWith("95"))) {
			return new String[] {famIndId.substring(0, 5), Integer.parseInt(famIndId.substring(5))+""};
		}

		if (famIndId.startsWith("ND")) {
			try {
				Integer.parseInt(famIndId.substring(2));
				return new String[] {famIndId, famIndId};
			} catch (NumberFormatException e) {}
		}

		if (index==-1) {
			if (famIndId.startsWith("7")) {
				System.err.println("Error - couldn't split '"+famIndId+"' into a FamilyID and an IndID");
			}
			return new String[] {famIndId, famIndId};
		}
		if (famIndId.toLowerCase().endsWith("a")||famIndId.toLowerCase().endsWith("b")) {
			System.err.println(famIndId);
			famIndId = famIndId.substring(0, famIndId.length()-1);
		}
		return new String[] {famIndId.substring(0, index), famIndId.substring(index+1)};
	}

	public static Hashtable<String,String> getBestPDdx() {
		BufferedReader reader = null;
		String[] line;
		String trav, prev;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = tools.getNinfoReader(2, false);

			ext.checkHeader(reader.readLine().split("[\\s]+"), NINFO2_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line.length>6) {
					if (ext.indexOfStr(line[6], PD_NAMES)>=0) {
						if (ext.indexOfStr(line[7], new String[] {"VPD", "CONF"})>=0) {
							trav = "VPD";
						} else if (ext.indexOfStr(line[7], new String[] {"NVPD", "POSS", "PROB"})>=0) {
							trav = "NVPD";
						} else if (line[7].equals("RPD")||line[7].equals("FRPT")||line[7].equals("NXAM")) {
							trav = "RPD";
						} else if (ext.indexOfStr(line[7], UNAFFECTED)>=0) {
							trav = line[7];
						} else {
							System.err.println("Error - Unknown PD designation for individual "+line[0]+"-"+line[1]+": '"+line[7]+"'");
							trav = "0";
						}
						prev = hash.containsKey(line[0]+"\t"+line[1])?hash.get(line[0]+"\t"+line[1]):"0";
						if (prev.equals("VPD")) {
							trav = "VPD";
						} else if (prev.equals("NVPD")) {
							trav = "NVPD";
						}
						hash.put(line[0]+"\t"+line[1], trav);
						hash.put(line[0]+ext.formNum(line[1], 3), trav);
					}
				} else {
					System.err.println("'"+line[0]+"-"+line[1]+"' shorted me");
				}
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error parsing ninfo2 file");
			System.exit(3);
		}

		return hash;
	}

	public static Hashtable<String,String> getExactPDdx() {
		return getExactPDdx(2);
	}

	public static Hashtable<String,String> getExactPDdx(int whichNinfo) {
		BufferedReader reader = null;
		String[] line;
		String trav, prev;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = tools.getNinfoReader(whichNinfo, false);

			ext.checkHeader(reader.readLine().split("[\\s]+"), NINFO2_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line.length>6) {
					if (ext.indexOfStr(line[6], PD_NAMES)>=0) {
						if (ext.indexOfStr(line[7], new String[] {"CONF"})>=0) {
							trav = "CONF_PD";
						} else if (ext.indexOfStr(line[7], new String[] {"VPD"})>=0) {
							trav = "VPD";
						} else if (ext.indexOfStr(line[7], new String[] {"NVPD", "POSS", "PROB"})>=0) {
							trav = "NVPD";
						} else if (line[7].equals("RPD")||line[7].equals("FRPT")||line[7].equals("NXAM")) {
							trav = "RPD";
						} else if (ext.indexOfStr(line[7], UNAFFECTED)>=0) {
							trav = line[7];
						} else {
							System.err.println("Error - Unknown PD designation for individual "+line[0]+"-"+line[1]+": '"+line[7]+"'");
							trav = "0";
						}
						prev = hash.containsKey(line[0]+"\t"+line[1])?hash.get(line[0]+"\t"+line[1]):"0";
						if (prev.equals("CONF_PD")) {
							trav = "CONF_PD";
						} else if (prev.equals("VPD")&&!trav.equals("CONF_PD")) {
							trav = "VPD";
						} else if (prev.equals("NVPD")&&!trav.equals("CONF_PD")) {
							trav = "NVPD";
						}
						hash.put(line[0]+"\t"+line[1], trav);
						hash.put(line[0]+ext.formNum(line[1], 3), trav);
					}
				} else {
					System.err.println("'"+line[0]+"-"+line[1]+"' shorted me");
				}
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error parsing ninfo"+whichNinfo+" file");
			System.exit(3);
		}

		try {
			reader = tools.getNinfoReader(3);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[0].toUpperCase().startsWith("PHENO")) {
					hash.put(line[1]+"\t"+line[2], line[3]);
					hash.put(line[1]+ext.formNum(line[2], 3), line[3]);
				}
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error reading file ninfo3");
			System.exit(2);
		}

		return hash;
	}

	public static boolean isAffected(Hashtable<String, String> hash, String trav) {
		String aff;

		// trav = trav.replace('\t',':');

		if (hash.containsKey(trav)) {
			aff = hash.get(trav);
			if (aff != null && (aff.equals("VPD")||aff.equals("NVPD")||aff.equals("RPD"))) {
				return true;
			}
		}

		return false;
	}

	public static String isAffected(String dx) {
		if (dx == null) {
			return ".";
		}
		if (dx.equals("CONF_PD")||dx.equals("VPD")||dx.equals("NVPD")||dx.equals("RPD")) {
			return "1";
		} else if (dx.equals("NOEV")||dx.equals("NRPD")) {
			return "0";
		} else {
			return ".";
		}
	}

	public static String isAffectedAndSeen(String dx) {
		if (dx == null) {
			return ".";
		}
		if (dx.equals("CONF_PD")||dx.equals("VPD")||dx.equals("NVPD")) {
			return "1";
		} else if (dx.equals("NOEV")) {
			return "0";
		} else {
			return ".";
		}
	}

	public static String isVPD(String dx) {
		if (dx == null) {
			return ".";
		}
		if (dx.equals("CONF_PD")||dx.equals("VPD")) {
			return "1";
		} else if (dx.equals("NVPD")) {
			return "0";
		} else {
			return ".";
		}
	}

	public static String isConfPD(String dx) {
		if (dx == null) {
			return ".";
		}
		if (dx.equals("CONF_PD")) {
			return "1";
		} else if (dx.equals(".")) {
			return ".";
		} else {
			return "0";
		}
	}

	public static Hashtable<String,String> getGender() {
		BufferedReader reader = null;
		String[] line;
		String trav;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = tools.getNinfoReader(2, false);

			ext.checkHeader(reader.readLine().split("[\\s]+"), NINFO2_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line[2].equals("F")) {
					trav = "2";
				} else if (line[2].equals("M")) {
					trav = "1";
				} else if (line[2].equals("N") || line[2].equals("?")) {
					trav = "0";
				} else {
					System.err.println("Error - unknown gender code for "+line[0]+"-"+line[1]+" ("+line[2]+")");
					trav = "0";
				}
				hash.put(line[0]+"\t"+line[1], trav);
			}
			reader.close();
		} catch (IOException ioe) {
			System.err.println("Error parsing ninfo2 file");
			System.exit(3);
		}

		return hash;
	}

	public static BufferedReader getNinfoReader(int whichNinfo) {
		return getNinfoReader(whichNinfo, true);
	}

	public static BufferedReader getNinfoReader(int whichNinfo, boolean kill) {
		String[] alt_locs = new String[] {tools.NINFO_DIR, ALT_DIR};
		BufferedReader reader = null;

		if (!ext.containsAny(whichNinfo+"", new String[] {"1", "2", "3", "5", "6"})) {
			System.err.println("Error - getNinfoReader only works for ninfos 1-3,6");
			System.exit(1);
		}

		try {
			reader = Files.getReader("ninfo"+whichNinfo+".dat", alt_locs);
			if (whichNinfo==1) {
				reader.mark(10000);
				ext.checkHeader(reader.readLine().trim().split("[\\s]+"), NINFO1_HEADER, true);
				reader.reset();
			}
			if (whichNinfo==2) {
				reader.mark(10000);
				ext.checkHeader(reader.readLine().trim().split("[\\s]+"), NINFO2_HEADER, true);
				reader.reset();
			}
		} catch (Exception e) {
			System.err.println("Error - ninfo"+whichNinfo+".dat not found in the current, 00ninfos, or 00masters directories");
			if (whichNinfo==3) {
				System.err.println("        --> assuming there are no autopsies, MZ twins, half sibs or ID-DNA# mixups that need to be addressed");
			} else {
				System.err.println("Please rectify");
			}
			if (kill) {
				System.exit(1);
			}
		}

		return reader;
	}

	public static Hashtable<String,String> pullTraitFromDB(String trait) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int index;
		String dbFile = DB_FILE;

		if (new File(dbFile).exists()) {
			System.out.println("Using "+dbFile+" in current directory");
		} else if (new File(CRF_DIR+dbFile).exists()) {
			dbFile = CRF_DIR+dbFile;
		} else if (new File(ALT_DIR+dbFile).exists()) {
			dbFile = ALT_DIR+dbFile;
		} else {
			System.err.println("Error - could not find the database ("+DB_FILE+") in current directory or alternate:\n"+ALT_DIR);
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(dbFile));
			index = ext.indexFactors(new String [] {trait}, reader.readLine().split("[\\s]+"), true, true)[0];
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				hash.put(line[1]+"\t"+line[2], line[index]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dbFile+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dbFile+"\"");
			System.exit(2);
		}

		return hash;
	}

	/**
	 * Returns the members of the largest sibship for each family
	 * 
	 * @return a matrix with a row for each family whith famid, father, mother, and the IndID of each member of the sibship
	 */
	public static String[][] getCombonentsOfTheMajorSibships(boolean onlyIncludeAffecteds, boolean onlyCountVPD) {
		BufferedReader reader;
        String[] line, fams, rents;
        String trav;
        Hashtable<String,Hashtable<String,Vector<String>>> famHash;
        Hashtable<String,Vector<String>> hash;
        Hashtable<String,String> bestDx;
        Vector<String> v;
        int count, max;
        String[][] results;
        int pick;
        
        bestDx = getBestPDdx();

        famHash = new Hashtable<String,Hashtable<String,Vector<String>>>();
        try {
	        reader = getNinfoReader(2, true);
			ext.checkHeader(reader.readLine().split("[\\s]+"), NINFO2_HEADER, true);
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (famHash.containsKey(line[0])) {
	        		hash = famHash.get(line[0]);
	        	} else {
	        		famHash.put(line[0], hash = new Hashtable<String,Vector<String>>());
	        	}
    			if (!onlyIncludeAffecteds || (onlyCountVPD && isVPD(bestDx.get(line[0]+"\t"+line[1])).equals("1")) || (!onlyCountVPD && isAffectedAndSeen(bestDx.get(line[0]+"\t"+line[1])).equals("1"))) {
    				HashVec.addToHashVec(hash, line[4]+"\t"+line[5], line[1], true);
    			}
	        }
	        reader.close();
        } catch (IOException ioe) {
	        System.err.println("Error parsing ninfo2.dat");
	        System.exit(2);
        }
        
        fams = HashVec.getKeys(famHash);
        results = new String[fams.length][];
        for (int i = 0; i<fams.length; i++) {
        	hash = famHash.get(fams[i]);
        	rents = HashVec.getKeys(hash);
        	pick = -1;
        	max = 0;
        	for (int j = 0; j<rents.length; j++) {
        		v = hash.get(rents[j]);
        		count = 0;
        		for (int k = 0; k<v.size(); k++) {
        			trav = v.elementAt(k);
        			if ((onlyCountVPD && isVPD(bestDx.get(fams[i]+"\t"+trav)).equals("1")) || (!onlyCountVPD && isAffectedAndSeen(bestDx.get(fams[i]+"\t"+trav)).equals("1"))) {
        				count++;
        			}
                }
        		if (count > max && !rents[j].split("[\\s]+")[0].equals("0") && !rents[j].split("[\\s]+")[1].equals("0")) {
        			max = count;
        			pick = j;
        		}
            }
        	if (pick == -1) {
        		results[i] = new String[] {fams[i]};
        	} else {
        		v = hash.get(rents[pick]);
        		results[i] = new String[3+v.size()];
        		results[i][0] = fams[i];
        		results[i][1] = rents[pick].split("[\\s]+")[0];
        		results[i][2] = rents[pick].split("[\\s]+")[1];
        		for (int j = 0; j < v.size(); j++) {
        			results[i][3+j] = v.elementAt(j);
				}
        	}
        }
        
        return results;
	}
}
