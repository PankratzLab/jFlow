package filesys;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

import common.*;

public class FamilyStructure {
	public static final String[][] TYPICAL_HEADERS = {{"FID", "famid"}, {"IID", "id"}, {"fa"}, {"mo"}, {"sex"}};

    private static final byte MISSING_VALUE_BYTE = (byte) -9;
    
    public static final int FID_INDEX = 0; 
    public static final int IID_INDEX = 1; 
    public static final int FA_INDEX = 2; 
    public static final int MO_INDEX = 3;

    public static final Object MISSING_ID_STR = "0"; 
	
    public ArrayList<String[]> cached_poPairsIDs = null;
    public ArrayList<int[]> cached_poPairsCompleteOnly = null;
    public HashMap<String, ArrayList<String>> cached_parentToChildrenMap = null;
    public boolean cached_parentMapIsCompleteOnly = false;
    public ArrayList<String[]> cached_all_trios = null;
    public ArrayList<int[]> cached_complete_trios = null;
    public ArrayList<String[]> cached_sib_pairs = null;
    public HashMap<String,Integer> cached_fidiidToIndexMap = null;
    
    public void clearCache() {
        cached_poPairsIDs = null;
        cached_poPairsCompleteOnly = null;
        cached_parentToChildrenMap = null;
        cached_parentMapIsCompleteOnly = false;
        cached_all_trios = null;
        cached_complete_trios = null;
        cached_sib_pairs  = null;
        cached_fidiidToIndexMap = null;
    }
    
	protected String[][] ids;
	protected String[] iids;
	protected String[] fids;
	protected String[] fas;
	protected String[] mos;
	protected byte[] genders;
	protected byte[] affections;
	protected String[] dnas;
	protected String[] mzTwinIds;
	protected Logger log;
	
	public FamilyStructure(String[][] ids, byte[] genders, byte[] affections) {
		this(ids, genders, affections, null);
	}
	
	public FamilyStructure(String[][] ids, byte[] genders, byte[] affections, String[] dnas) {
	    this(ids, genders, affections, dnas, new String[ids.length], new Logger());
	}
	
	public FamilyStructure(String[][] ids, byte[] genders, byte[] affections, String[] dnas, String[] mzTwinIds, Logger log) {
		this.ids = ids;
		this.genders = genders;
		this.affections = affections;
		this.dnas = dnas;
		this.mzTwinIds = mzTwinIds;
		this.log = log;
	}
	
	public FamilyStructure(String filename) {
		this(filename, false);
	}
	
	public FamilyStructure(String filename, boolean loadDNAs) {
	    this(filename, loadDNAs, new Logger());
	}

	public FamilyStructure(String filename, boolean loadDNAs, Logger log) {
		BufferedReader reader;
		String[] line;
		int count;

		try {
			count = Files.countLines(filename, 0);
			reader = Files.getAppropriateReader(filename);
			ids = new String[count][];
			iids = new String[count];
			fids = new String[count];
			fas = new String[count];
			mos = new String[count];
			genders = new byte[count];
			affections = new byte[count];
			dnas = loadDNAs ? new String[count] : null;
			mzTwinIds = new String[count];
			for (int i = 0; i < count; i++) {
				line = reader.readLine().trim().split("[\\s]+");
				ids[i] = new String[] {line[0], line[1], line[2], line[3]};
				fids[i] = line[0];
				iids[i] = line[1];
				fas[i] = line[2];
				mos[i] = line[3];
				genders[i] = ext.isMissingValue(line[4]) ? FamilyStructure.MISSING_VALUE_BYTE : Byte.parseByte(line[4]);  
				affections[i] = ext.isMissingValue(line[5]) ? FamilyStructure.MISSING_VALUE_BYTE : Byte.parseByte(line[5]); 
				if (loadDNAs) {
					dnas[i] = line[6]; 
				}
				if(line.length > 7){
					mzTwinIds[i] = ext.isMissingValue(line[7]) ? null : line[7];
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
		    this.log.reportFileNotFound(filename);
		} catch (IOException ioe) {
		    this.log.reportIOException(filename);
		}
	}
	
	public String[][] getIDs() {
		return ids;
	}
	
	public int getIndIndex(String fid, String iid){
		Integer index = getfidiidToIndexMap().get(fid + "\t" + iid);
		return index == null ? -1 : index;
	}
	
	public HashMap<String,Integer> getfidiidToIndexMap() {
		if (cached_fidiidToIndexMap == null) {
			buildfidiidToIndexMap();
		}
		return cached_fidiidToIndexMap;
	}
	
	private void buildfidiidToIndexMap(){
		cached_fidiidToIndexMap = new HashMap<String,Integer>();
		for(int i = 0; i < ids.length; i++){
			if (cached_fidiidToIndexMap.put(this.fids[i] + "\t" + this.iids[i], i) != null){
				System.err.println("Warning - Pedigree contains non-unique FID/IID combinations!");
			}
		}
	}
	
    public String getiDNA(int i) {
        return this.dnas[i];
    }
	
	public String getFID(int index) {
	    return this.ids[index][FID_INDEX];
	}
	
	public String getIID(int index) {
	    return this.ids[index][IID_INDEX];
	}
	
	public String getFA(int index) {
	    return this.ids[index][FA_INDEX];
	}
	
	public int getIndexOfFaInIDs(int indivIndex) {
		return MISSING_ID_STR.equals(fas[indivIndex]) ? -1 : getIndIndex(fids[indivIndex],fas[indivIndex]);
	}
	
	public int getIndexOfMoInIDs(int indivIndex) {
		return MISSING_ID_STR.equals(mos[indivIndex]) ? -1 : getIndIndex(fids[indivIndex],mos[indivIndex]);
	}
	
	public String getMO(int index) {
	    return this.ids[index][MO_INDEX];
	}
	
	public byte[] getGenders() {
		return genders;
	}
	
	public byte getGender(int index) {
	    return genders[index];
	}

	public byte[] getAffections() {
		return affections;
	}

	public void setAffections(byte[] affs) {
		affections = affs;
	}
	
	public String[] getDnas() {
		return dnas;
	}
	
	public String[] getMzTwinIds() {
		return mzTwinIds;
	}
	
	public String getMzTwinId(int index) {
		return mzTwinIds[index];
	}
	
	public String getIndividualHeader(int index, boolean displayDNA) {
		return Array.toStr(ids[index]) + "\t" + genders[index] + "\t" + affections[index] + (displayDNA && dnas != null ? "\t" + dnas[index] : "");
	}

	public void writeToFile(String filename, boolean displayDNA) {
        PrintWriter writer;
        
        try {
	        writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i<ids.length; i++) {
				writer.println(getIndividualHeader(i, displayDNA));
	        }
	        writer.close();
        } catch (Exception e) {
            this.log.reportException(e);
        }
	}
	
	public static boolean likelyPedHeader(String[] line) {
		return Array.countIf(ext.indexFactors(TYPICAL_HEADERS, line, false, true, false, false), -1) < 3;
	}
}
