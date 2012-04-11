package cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Properties;
import java.util.Vector;

import common.Array;
import common.Files;
import common.HashVec;
import common.ext;
import cnv.var.SampleData;

public class Project extends Properties {
	public static final long serialVersionUID = 1L;
//	public static final String DEFAULT_PROJECT = "these.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/pd_win.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/load_win.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/boss.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/HumanHap550_win.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/consortium.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/consortiumReplication.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/consortiumReplicationLRR.properties";
	public static final String DEFAULT_PROJECT = "/Users/zxu/workspace/Genvisis/projects/practice.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/strat_demo.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/boss.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/load_win.properties";
//	public static final String DEFAULT_PROJECT = "/home/npankrat/projects/demo.proj";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/pd_win.properties";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/HumanHap550_win.properties";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/myron_excision.proj";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/demo_excision.proj";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/sing550_win.proj";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/consortium.properties";
	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/consortiumReplication.properties";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/demo_indian_diabetes.proj";
//	public static final String DEFAULT_CURRENT = "/home/npankrat/projects/boss.proj";
	
	public static final String DEFAULT_SCATTER_PROJECT = "/home/npankrat/projects/demo.properties";
//	public static final String DEFAULT_SCATTER_PROJECT = "/home/npankrat/projects/pd_win.properties";
//	public static final String DEFAULT_SCATTER_PROJECT = "/home/npankrat/projects/load_win.properties";
//	public static final String DEFAULT_SCATTER_PROJECT = "/home/npankrat/projects/sing550_win.proj";


	public static final String DEFAULT_PROPERTIES = "cnv/filesys/default.properties";

	public static final String PROJECT_NAME = "PROJECT_NAME";
	public static final String PROJECT_PROPERTIES_FILENAME = "FILENAME";
	public static final String PROJECT_DIRECTORY = "PROJECT_DIRECTORY";
	public static final String SOURCE_DIRECTORY = "SOURCE_DIRECTORY";
	public static final String SOURCE_FILENAME_EXTENSION = "SOURCE_FILENAME_EXTENSION";
	public static final String ID_HEADER = "ID_HEADER";
	public static final String PARSE_AT_AT_SYMBOL = "PARSE_AT_AT_SYMBOL";
	public static final String JAR_STATUS = "JAR_STATUS";
	public static final String SAMPLE_DIRECTORY = "SAMPLE_DIRECTORY";
	public static final String IND_DIRECTORY = "IND_DIRECTORY";
	public static final String DATA_DIRECTORY = "DATA_DIRECTORY";
	public static final String PLOT_DIRECTORY = "PLOT_DIRECTORY";
	public static final String RESULTS_DIRECTORY = "RESULTS_DIRECTORY";
	public static final String DEMO_DIRECTORY = "DEMO_DIRECTORY";
	public static final String MARKER_POSITION_FILENAME = "MARKER_POSITION_FILENAME";
	public static final String MARKERSET_FILENAME = "MARKERSET_FILENAME";
	public static final String MARKERLOOKUP_FILENAME = "MARKERLOOKUP_FILENAME";
	public static final String SAMPLELIST_FILENAME = "SAMPLELIST_FILENAME";
	public static final String SAMPLE_SUBSET_FILENAME = "SAMPLE_SUBSET_FILENAME";
	public static final String SAMPLE_DATA_FILENAME = "SAMPLE_DATA_FILENAME";
	public static final String ORIGINAL_CENTROIDS_FILENAME = "ORIGINAL_CENTROIDS_FILENAME";
	public static final String GENOTYPE_CENTROIDS_FILENAME = "GENOTYPE_CENTROIDS_FILENAME";
	public static final String CHIMERA_CENTROIDS_FILENAME = "CHIMERA_CENTROIDS_FILENAME";
	public static final String CUSTOM_CENTROIDS_FILENAME = "CUSTOM_CENTROIDS_FILENAME";
	public static final String DISPLAY_MARKERS_FILENAME = "DISPLAY_MARKERS_FILENAME";
	public static final String FILTERED_MARKERS_FILENAME = "FILTERED_MARKERS_FILENAME";
	public static final String PEDIGREE_FILENAME = "PEDIGREE_FILENAME";
	public static final String MOSAIC_COLOR_CODES_FILENAME = "MOSAIC_COLOR_CODES_FILENAME";
	public static final String MOSAIC_RESULTS_FILENAME = "MOSAIC_RESULTS_FILENAME";
	public static final String REGION_LIST_FILENAMES = "REGION_LIST_FILENAMES";
	public static final String CNV_FILENAMES = "CNV_FILENAMES";
	public static final String STRATIFICATION_RESULTS_FILENAMES = "STRATIFICATION_RESULTS_FILENAMES";
	public static final String GC_THRESHOLD = "GC_THRESHOLD";
	public static final String QQ_FILENAMES = "QQ_FILENAMES";
	public static final String DISPLAY_QUANTILES = "DISPLAY_QUANTILES";
	public static final String DISPLAY_STANDARD_QQ = "DISPLAY_STANDARD_QQ";
	public static final String DISPLAY_ROTATED_QQ = "DISPLAY_ROTATED_QQ";
	public static final String TARGET_MARKERS_FILENAME = "TARGET_MARKERS_FILENAME";
	public static final String SOURCE_FILE_DELIMITER = "SOURCE_FILE_DELIMITER";
	public static final String PENNCNV_EXECUTABLE_DIRECTORY = "PENNCNV_EXECUTABLE_DIRECTORY";
	public static final String PENNCNV_DATA_DIRECTORY = "PENNCNV_DATA_DIRECTORY";
	public static final String PENNCNV_RESULTS_DIRECTORY = "PENNCNV_RESULTS_DIRECTORY";
	public static final String LRRSD_CUTOFF = "LRRSD_CUTOFF";
	public static final String MOSAIC_ARMS_FILENAME = "MOSAIC_ARMS_FILENAME";
	public static final String NUM_THREADS = "NUM_THREADS";
	public static final String LONG_FORMAT = "LONG_FORMAT";
	public static final String NUM_MARKERS_PER_FILE = "NUM_MARKERS_PER_FILE";
	public static final String CLUSTER_FILTER_COLLECTION_FILENAME = "CLUSTER_FILTER_COLLECTION_FILENAME";
	public static final String SEXCHECK_RESULTS_FILENAME = "SEXCHECK_RESULTS_FILENAME";

	private boolean jar;

	public Project() {
		Files.loadProperties(this, DEFAULT_PROPERTIES, true, true, false);
	}
	
	public Project(String filename, boolean jar) {
		this();
		
		screenProperties(filename);
		Files.loadProperties(this, filename, jar, true, false);

        setProperty(PROJECT_DIRECTORY, ext.verifyDirFormat(getProperty(PROJECT_DIRECTORY)));
        setProperty(SOURCE_DIRECTORY, ext.verifyDirFormat(getProperty(SOURCE_DIRECTORY)));
        setProperty(PROJECT_PROPERTIES_FILENAME, filename);
        
        this.jar = jar;
	}

	public String getDir(String directory) {
		return getDir(directory, false);
	}
	
	public String getDir(String directory, boolean mkdirs) {
		String dir = null;
		
		if (containsKey(directory)) {
			dir = getProperty(directory);
			if (!dir.startsWith("/") && dir.indexOf(":") == -1) {
				dir = getProperty(PROJECT_DIRECTORY)+dir;
			}
			if (!Files.exists(dir, jar)) {
				if (mkdirs && !jar) {
					new File(dir).mkdirs();
				} else {
					System.err.println("Error - file '"+dir+"' does not exist");
				}
			}
		} else {
			System.err.println("Error - directory '"+directory+"' is undefined in cnv.filesys.Project");
		}
		
		return dir;
	}
	
	public String getProjectDir() {
		return getProperty(PROJECT_DIRECTORY);
	}

	public String getFilename(String fileType) {
		return getFilename(fileType, false, true);
	}
	
	public String getFilename(String fileType, boolean make, boolean verbose) {
		return getFilename(fileType, null, make, verbose);
	}
	
	public String getFilename(String fileType, String subdirectory, boolean make, boolean verbose) {
		String file = null;
		
		if (containsKey(fileType)) {
			file = getProperty(fileType);
			if (!file.startsWith("/") && file.indexOf(":") == -1) {
				file = getProperty(PROJECT_DIRECTORY)+(subdirectory==null?"":getProperty(subdirectory))+file;
			}
			if (!Files.exists(file, getJarStatus())) {
				if (make) {
					new File(ext.parseDirectoryOfFile(file)).mkdirs();
				} else if (verbose) {
					System.err.println("Error - file '"+file+"' does not exist");
				}
			}
		} else {
			System.err.println("Error - file '"+fileType+"' is undefined in cnv.filesys.Project");
		}

		return file;
	}
	
	public String[] getFilenames(String type) {
		return getFilenames(type, false);
	}
	
	public String[] getFilenames(String type, boolean suppressMissing) {
		String[] files = null;
		Vector<String> v;
		
		v = new Vector<String>();
		if (containsKey(type)) {
			files = getProperty(type).split(";");
			if (files.length == 1 && files[0].equals("")) {
				files = new String[0];
			}
			
			for (int i = 0; i<files.length; i++) {
				if (!files[i].startsWith("/") && files[i].indexOf(":") == -1) {
					files[i] = getProperty(PROJECT_DIRECTORY)+files[i];
				}
				if (!Files.exists(files[i], getJarStatus()) && !suppressMissing) {
					System.err.println("Error - file '"+files[i]+"' does not exist");
				} else {
					v.add(files[i]);
				}
            }
		} else {
			System.err.println("Error - file '"+type+"' is undefined in cnv.filesys.Project");
		} 

		return Array.toStringArray(v);
	}
	
	public double getDouble(String variable) {
		String trav;
		
		trav = getProperty(variable);
		try {
			return Double.parseDouble(trav);
		} catch (NumberFormatException nfe) {
			System.err.println("Error - '"+trav+"' is not a valid value for "+variable);
			return Double.NaN;
		}		
	}
	

	public MarkerSet getMarkerSet() {
		if (Files.exists(getFilename(MARKERSET_FILENAME), getJarStatus())) {
			return MarkerSet.load(getFilename(MARKERSET_FILENAME), getJarStatus());
		} else {
			return null;
		}
	}

	public String[] getMarkerNames() {
		return getMarkerSet().getMarkerNames();
	}

	public MarkerLookup getMarkerLookup() {
		if (Files.exists(getFilename(MARKERLOOKUP_FILENAME), getJarStatus())) {
			return MarkerLookup.load(getFilename(MARKERLOOKUP_FILENAME), getJarStatus());
		} else {
			return null;
		}
	}

	public SampleList getSampleList() {
		if (Files.exists(getFilename(SAMPLELIST_FILENAME), getJarStatus())) {
			return SampleList.load(getFilename(SAMPLELIST_FILENAME), getJarStatus());
		} else {
			System.out.println("Failed to find SampleList; generating one...");
			return SampleList.generateSampleList(this);
		}
	}

	public String[] getSamples() {
		SampleList sampleList;
		
		sampleList = getSampleList();
		if (sampleList == null) {
			return null;
		} else {
			return sampleList.getSamples();
		}
	}

	public FullSample getFullSample(String sample) {
		if (Files.exists(getDir(SAMPLE_DIRECTORY)+sample+".fsamp", getJarStatus())) {
			return FullSample.load(getDir(SAMPLE_DIRECTORY)+sample+".fsamp", getJarStatus());
		} else {
			return null;
		}
	}

	public Sample getSample(String sample) {
		if (Files.exists(getDir(IND_DIRECTORY)+sample+".samp", getJarStatus())) {
			return Sample.load(getDir(IND_DIRECTORY)+sample+".samp", getJarStatus());
		} else {
			return null;
		}
	}

	public SampleData getSampleData(boolean loadCNVs) {
		return new SampleData(this, loadCNVs);
	}

	public Hashtable<String,String> getFilteredHash() {
		if (getProperty(FILTERED_MARKERS_FILENAME).equals("")) {
			return new Hashtable<String,String>();
		} else if (Files.exists(getFilename(FILTERED_MARKERS_FILENAME), getJarStatus())) {
			return HashVec.loadFileToHashString(getFilename(FILTERED_MARKERS_FILENAME), 0, new int[] {0}, "", false, getJarStatus());
		} else {
			System.err.println("Error - '"+getProperty(FILTERED_MARKERS_FILENAME)+"' not found");
			return new Hashtable<String,String>();
		}
	}


	public Vector<String> getStratResults() {
		String[] files;
		Vector<String> v;

		files = Files.list(getProjectDir(), ".mds", getJarStatus());
	
		v = new Vector<String>();
		if (files == null) {
			System.err.println("Error - no .mds files found in directory");
		} else {
			for (int i = 0; i<files.length; i++) {
				v.add(getProjectDir()+files[i]);
				System.out.println(getProjectDir()+files[i]);
	        }
		}
		
		return v;
	}
	
	public String getSourceFileDelimiter() {
		String str;
		
		str = getProperty(SOURCE_FILE_DELIMITER);
		
		if (str.toUpperCase().equals("COMMA")) {
			return ",";
		} else if (str.toUpperCase().equals("TAB")) {
			return "\t";
		} else if (str.toUpperCase().equals("SPACE")) {
			return "[\\s]+";
		} else {
			System.err.println("Error - invalid delimiter specified: '"+str+"'");
			return ",";
		} 
	}	
	
	public boolean getJarStatus() {
		return jar;
	}
	
	public String[] getRegionLists() {
		return getProperty(REGION_LIST_FILENAMES).trim().split(";");
	}
	
	public void screenProperties(String filename) {
		BufferedReader reader;
        String[] line;
        String trav;
        boolean changed;
        Vector<String> knowns, unknowns, corrections;
        
        changed = false;
        knowns = Array.toStringVector(HashVec.getKeys(this, false, false));
        unknowns = new Vector<String>();
        corrections = new Vector<String>();
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        while (reader.ready()) {
	        	trav = reader.readLine();
	        	line = trav.trim().split("=");
	        	if (line[0].startsWith("#") || line[0].startsWith("??_") || line[0].equals("")) {
	        		corrections.add(trav);
	        		if (line[0].length() > 0 && getProperty(line[0].substring(1)) != null) {
	        			knowns.remove(line[0].substring(1));
	        		}
	        	} else {
	        		if (getProperty(line[0]) == null) {
	        			unknowns.add(line[0]);
	        			corrections.add("??_"+trav);
	        		} else {
	        			knowns.remove(line[0]);
	        			corrections.add(trav);
	        		}
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
        
        if (unknowns.size() > 0) {
        	System.err.println("Error - check spelling for the following unexpected propert"+(unknowns.size()==1?"y":"ies")+" in "+filename+":");
        	changed = true;
        }
        for (int i = 0; i < unknowns.size(); i++) {
        	System.err.println("        "+unknowns.elementAt(i));
		}
        
        if (knowns.size() > 0) {
        	changed = true;
        	corrections.add("");
        	corrections.add("# A few more paramters that were not originally defined:");
        }
        for (int i = 0; i < knowns.size(); i++) {
        	corrections.add("#"+knowns.elementAt(i)+"="+getProperty(knowns.elementAt(i)));
		}        
        
        if (changed) {
        	new File(filename).renameTo(new File(filename+".bak"));
        	Files.backup(ext.removeDirectoryInfo(filename+".bak"), ext.parseDirectoryOfFile(filename), ext.parseDirectoryOfFile(filename), true);
        	Files.writeList(Array.toStringArray(corrections), filename);
        }
	}

	public static void main(String[] args) {

	}
}
