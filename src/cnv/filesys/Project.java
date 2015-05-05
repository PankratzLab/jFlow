package cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.manage.TransposeData;
import cnv.var.SampleData;
import filesys.GeneSet;

public class Project {

	private static String[] verifyFiles(String[] strs, boolean isDir) {
		String[] verified = new String[strs.length];
		for (int i = 0; i < strs.length; i++) {
			verified[i] = isDir ? ext.verifyDirFormat(strs[i]) : ext.replaceAllWith(strs[i], "\\", "/");
		}
		return verified;
	}
	
	public abstract static class Property<T> {
		private final Project myProj;
		private final String name;
		private final String desc;
		private final T defaultValue;
		private T value;
		public Property(Project proj, String name, String description, T defVal) {
			this.myProj = proj;
			this.name = name;
			this.desc = description;
			this.defaultValue = defVal;
			this.value = defaultValue;
		}
		public Project      getProject() { return myProj; }
		public       T        getValue() { return value; }
		public    void setValue(T value) { this.value = value; }
		public  String         getName() { return name; }
		public  String  getDescription() { return desc; }
		public       T getDefaultValue() { return defaultValue; }
		public  abstract void parseValue(String valueStr);
		public  String getValueString() { return value.toString(); }
		public  String getDefaultValueString() { return defaultValue.toString(); }
	}
	
	public static class StringProperty extends Property<String> {
		public StringProperty(Project proj, String name, String desc, String defVal) {
			super(proj, name, desc, defVal);
		}
		@Override
		public void parseValue(String valueStr) { setValue(valueStr); }
	}
	
	public static class StringListProperty extends Property<String[]> {
		static String delim = ";";
		public StringListProperty(Project proj, String name, String desc, String[] defVal) {
			super(proj, name, desc, defVal);
		}
		public StringListProperty(Project proj, String name, String desc, String defVal) {
			super(proj, name, desc, defVal.split(delim));
		}
		@Override
		public void parseValue(String valueStr) {
			this.setValue(valueStr.split(delim));
		}
		@Override
		public String getValueString() {
			return Array.toStr(getValue(), ";");
		}
		public String getDefaultValueString() { return Array.toStr(getDefaultValue(), ";"); }
	}
	
	public static class BooleanProperty extends Property<Boolean> {
		public BooleanProperty(Project proj, String name, String desc, Boolean defVal) {
			super(proj, name, desc, defVal);
		}
		@Override
		public void parseValue(String valueStr) { 
			Boolean newValue = Boolean.valueOf(valueStr);
			setValue(newValue); 
		}
	}
	
	public static class IntegerProperty extends Property<Integer> {
		int currValue;
		int min, max;
		public IntegerProperty(Project proj, String name, String desc, int min, int max, int defValue) {
			super(proj, name, desc, defValue);
			if (min > max || defValue < min || defValue > max || (max == min && defValue != max)) {
				throw new RuntimeException("Cannot initialize IntegerProperty with: min=" + min + ", max=" + max + ", and default value=" + defValue);
			}
			this.min = min;
			this.max = max;
		}
		@Override
		public void parseValue(String valueStr) {
			Integer newValue = Integer.valueOf(valueStr);
			setValue(newValue);
		}
		@Override
		public void setValue(Integer value) {
			if (value < min || value > max) {
				throw new RuntimeException("Error - values for property " + getName() + " must be within " + min + "-" + max + "; " + value + " is not valid");
			}
			super.setValue(value);
		}
	}
	
	public static class DoubleProperty extends Property<Double> {
		double currValue;
		double min;
		double max;
		public DoubleProperty(Project proj, String name, String desc, double min, double max, double defValue) {
			super(proj, name, desc, defValue);
			if (min > max || defValue < min || defValue > max || (max == min && defValue != max)) {
				throw new RuntimeException("Cannot initialize DoubleProperty['" + name + "'] with: min=" + min + ", max=" + max + ", and default value=" + defValue);
			}
			this.min = min;
			this.max = max;
		}
		@Override
		public void parseValue(String valueStr) {
			Double newValue = Double.valueOf(valueStr);
			setValue(newValue);
		}
		@Override
		public void setValue(Double value) {
			if (value < min || value > max) {
				throw new RuntimeException("Error - values for property " + getName() + " must be within " + min + "-" + max + "; " + value + " is not valid");
			}
			super.setValue(value);
		}
	}
	
	public static class FileProperty extends StringProperty {
		final boolean isDir;
		public FileProperty(Project proj, String name, String desc, String defVal, boolean dirOnly) {
			super(proj, name, desc, dirOnly ? ext.verifyDirFormat(defVal) : ext.replaceAllWith(defVal, "\\", "/")/* == null || "".equals(defVal) ? null : new File(defVal)*/);
			isDir = dirOnly;
		}
		@Override
		public void setValue(String value) {
			super.setValue(isDir ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/"));
		}
		@Override
		public String getValue() {
			return getValue(false, false);
		}
		public String getValue(boolean mkdirs, boolean verbose) {
			return getValue(null, mkdirs, verbose);
		}
		public String getValue(String subdir, boolean mkdirs, boolean verbose) {
			String valu = super.getValue();
			if (valu.contains("~")) {
				valu = ext.replaceTilde(valu);
			}
			
			String tempValue = valu;
			if (!"".equals(valu) /*&& !valu.startsWith(".")*/ && !valu.startsWith("/") && valu.indexOf(":") == -1) {
				if (isDir) {
					tempValue = getProject().PROJECT_DIRECTORY.getValue() + valu;
				} else {
					tempValue = getProject().PROJECT_DIRECTORY.getValue() + (subdir == null ? "" : getProject().getProperty(subdir).getValueString()) + valu;
				}
			}
			if (!Files.exists(tempValue, getProject().JAR_STATUS.getValue())) {
				if (mkdirs/* && getProject().JAR_STATUS.getValue()*/) {
					if (isDir) {
						(new File(tempValue)).mkdirs();
					} else {
						(new File(ext.parseDirectoryOfFile(tempValue))).mkdirs();
					}
				} else if (verbose) {
					if (isDir) {
						getProject().getLog().reportError("Error - directory '"+valu+"' does not exist");
					} else {
						if (!Files.exists(ext.parseDirectoryOfFile(valu), getProject().JAR_STATUS.getValue())) {
							getProject().getLog().reportError("Error - the directory ('"+ext.parseDirectoryOfFile(valu)+"') of the file you're trying to access/create ('"+ext.removeDirectoryInfo(valu)+"') does not exist");
						} else {
							getProject().getLog().reportError("Error - file '"+valu+"' does not exist");
						}
					}
				}
			}
			
			return tempValue;
		}
	}
	
	public static class MultiFileProperty extends Property<String[]> {
		final boolean isDir;
		public MultiFileProperty(Project proj, String name, String desc, String[] defVal, boolean dirOnly) {
			super(proj, name, desc, verifyFiles(defVal, dirOnly)/*createFiles(defVal)*/);
			isDir = dirOnly;
		}
		public MultiFileProperty(Project proj, String name, String desc, String defVal, boolean dirOnly) {
			super(proj, name, desc, verifyFiles(defVal.split(";"), dirOnly) /*createFiles(defVal)*/);
			isDir = dirOnly;
		}
		@Override
		public void setValue(String[] value) {
			super.setValue(verifyFiles(value, isDir));
		}
		@Override
		public void parseValue(String valueStr) {
			String[] newValue = valueStr.split(";");
			setValue(newValue);
		}
		@Override
		public String getValueString() {
			return Array.toStr(getValue(), ";");
		}
		public String getDefaultValueString() { return Array.toStr(getDefaultValue(), ";"); }

		@Override
		public String[] getValue() {
			return getValue(true);
		}
		
		public String[] getValue(boolean suppressMissing) {
			String[] files = super.getValue();
			Vector<String> v = new Vector<String>();

			if (files.length == 1 && files[0].equals("")) {
				files = new String[0];
			}
			
			for (int i = 0; i < files.length; i++) {
				if (null != files[i] && !"".equals(files[i]) /*&& !files[i].startsWith(".")*/ && !files[i].startsWith("/") && files[i].indexOf(":") == -1) {
					files[i] = getProject().PROJECT_DIRECTORY.getValue() + (files[i].startsWith("./") ? files[i].substring(2) : files[i]);
				}
				if (!Files.exists(files[i], getProject().JAR_STATUS.getValue()) && !suppressMissing) {
					System.err.println("Error - file '"+files[i]+"' does not exist");
				} //else {
				v.add(files[i]);
				//}
            }
			
			return Array.toStringArray(v);
		}
		
		
	}
	
	public static class EnumProperty extends StringProperty {
		String[] options;
		public EnumProperty(Project proj, String name, String desc, int defaultIndex, String... opts) {
			this(proj, name, desc, opts[defaultIndex], opts);
		}
		public EnumProperty(Project proj, String name, String desc, String defVal, String... opts) {
			super(proj, name, desc, defVal);
			this.options = opts;
			boolean found = false;
			for (String opt : this.options) {
				if (opt.equals(defVal)) {
					found = true;
					break;
				}
			}
			if (!found) {
				throw new RuntimeException("Cannot initialize EnumProperty['" + name + "'] with: default value=" + defVal + ", and options [" + Arrays.toString(opts) + "]");
			}
		}
		public String[] getValidOptions() { return options; }
	}
	
	public <T> T getProperty(Property<T> prop) {
		return prop.getValue();
	}
	
	public <T> void setProperty(Property<T> prop, T value) {
		prop.setValue(value);
	}
	
	public void setProperty(String name, String value) {
		getProperty(name).parseValue(value);
	}
	
	@SuppressWarnings("unchecked")
	public <T extends Property<?>> T getProperty(String name) {
		try {
			return (T) this.getClass().getField(name).get(this);
		} catch (Exception e) {
			return null;
//			throw new RuntimeException(e);
		}
	}
	
	public boolean containsKey(String name) {
		try {
			return this.getClass().getField(name).get(this) instanceof Property;
		} catch (SecurityException e) {
			throw new RuntimeException(e);
		} catch (Exception e) {
			return false;
		} 
	}
	
	public   IntegerProperty                            LOG_LEVEL = new   IntegerProperty(this,                            "LOG_LEVEL", "", -1, 12, 1);
	public    StringProperty                         PROJECT_NAME = new    StringProperty(this,                         "PROJECT_NAME", "Project Name", "New Project");
	public    StringProperty            SOURCE_FILENAME_EXTENSION = new    StringProperty(this,            "SOURCE_FILENAME_EXTENSION", "", ".csv");
	public    StringProperty                            ID_HEADER = new    StringProperty(this,                            "ID_HEADER", "", "Sample Name");
	public    StringProperty                TWOD_LOADED_VARIABLES = new    StringProperty(this,                "TWOD_LOADED_VARIABLES", "", "");
	public    StringProperty                            FID_ALIAS = new    StringProperty(this,                            "FID_ALIAS", "", "FID;F_ID;FamID;Fam_ID;Family;FamilyID;Family_ID");
	public    StringProperty                            IID_ALIAS = new    StringProperty(this,                            "IID_ALIAS", "", "ID;IID;I_ID;IndID;Ind_ID");
	public    StringProperty                         SAMPLE_ALIAS = new    StringProperty(this,                         "SAMPLE_ALIAS", "", "Sample;DNA;DNA#");
	public   BooleanProperty                   PARSE_AT_AT_SYMBOL = new   BooleanProperty(this,                   "PARSE_AT_AT_SYMBOL", "", Boolean.FALSE);
	public   BooleanProperty                           JAR_STATUS = new   BooleanProperty(this,                           "JAR_STATUS", "", Boolean.FALSE);
	public   BooleanProperty                    DISPLAY_QUANTILES = new   BooleanProperty(this,                    "DISPLAY_QUANTILES", "", Boolean.FALSE);
	public   BooleanProperty                  DISPLAY_STANDARD_QQ = new   BooleanProperty(this,                  "DISPLAY_STANDARD_QQ", "", Boolean.TRUE);
	public   BooleanProperty                   DISPLAY_ROTATED_QQ = new   BooleanProperty(this,                   "DISPLAY_ROTATED_QQ", "", Boolean.FALSE);
	public   BooleanProperty                   PENNCNV_GZIP_YESNO = new   BooleanProperty(this,                   "PENNCNV_GZIP_YESNO", "", Boolean.TRUE);
	public   BooleanProperty                          LONG_FORMAT = new   BooleanProperty(this,                          "LONG_FORMAT", "", Boolean.FALSE);
	public   BooleanProperty           SHIFT_SEX_CHR_COLORS_YESNO = new   BooleanProperty(this,           "SHIFT_SEX_CHR_COLORS_YESNO", "", Boolean.TRUE);
	public    DoubleProperty                         GC_THRESHOLD = new    DoubleProperty(this,                         "GC_THRESHOLD", "", 0.0, 1.0, 0.15);
	public    DoubleProperty                         LRRSD_CUTOFF = new    DoubleProperty(this,                         "LRRSD_CUTOFF", "", 0.0, 1.0, 0.32);
	public   IntegerProperty                          NUM_THREADS = new   IntegerProperty(this,                          "NUM_THREADS", "", 1, 99, 1);
	public   IntegerProperty              QQ_MAX_NEG_LOG10_PVALUE = new   IntegerProperty(this,              "QQ_MAX_NEG_LOG10_PVALUE", "", 1, 10000, 100);
	public   IntegerProperty WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER = new   IntegerProperty(this, "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER", "", 1, 1000000, 10000);
	public   IntegerProperty         MAX_MARKERS_LOADED_PER_CYCLE = new   IntegerProperty(this,         "MAX_MARKERS_LOADED_PER_CYCLE", "", 1, 10000, 100);
	public   IntegerProperty  MAX_MEMORY_USED_TO_LOAD_MARKER_DATA = new   IntegerProperty(this,  "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA", "", 8, 65536, 250);
	public   IntegerProperty          INTENSITY_PC_NUM_COMPONENTS = new   IntegerProperty(this,          "INTENSITY_PC_NUM_COMPONENTS", "", 0, 10000, 100);
	public      FileProperty                    PROJECT_DIRECTORY = new      FileProperty(this,                    "PROJECT_DIRECTORY", "", "./", true);
	public      FileProperty                     SOURCE_DIRECTORY = new      FileProperty(this,                     "SOURCE_DIRECTORY", "", "./", true);
	public      FileProperty                     SAMPLE_DIRECTORY = new      FileProperty(this,                     "SAMPLE_DIRECTORY", "", "./samples/", true);
	public      FileProperty                       DATA_DIRECTORY = new      FileProperty(this,                       "DATA_DIRECTORY", "", "./data/", true);
	public      FileProperty                MARKER_DATA_DIRECTORY = new      FileProperty(this,                "MARKER_DATA_DIRECTORY", "", "./transposed/", true);
	public      FileProperty                    RESULTS_DIRECTORY = new      FileProperty(this,                    "RESULTS_DIRECTORY", "", "./results/", true);
	public      FileProperty                       DEMO_DIRECTORY = new      FileProperty(this,                       "DEMO_DIRECTORY", "", "./demo/", true);
	public      FileProperty         PENNCNV_EXECUTABLE_DIRECTORY = new      FileProperty(this,         "PENNCNV_EXECUTABLE_DIRECTORY", "", "/home/npankrat/bin/", true);
	public      FileProperty               PENNCNV_DATA_DIRECTORY = new      FileProperty(this,               "PENNCNV_DATA_DIRECTORY", "", "./penn_data/", true);
	public      FileProperty            PENNCNV_RESULTS_DIRECTORY = new      FileProperty(this,            "PENNCNV_RESULTS_DIRECTORY", "", "./penncnv/", true);
	public      FileProperty                     BACKUP_DIRECTORY = new      FileProperty(this,                     "BACKUP_DIRECTORY", "", "./backup/", true);
	public      FileProperty          PROJECT_PROPERTIES_FILENAME = new      FileProperty(this,                             "FILENAME", "", "./example.properties", false);
	public      FileProperty             MARKER_POSITION_FILENAME = new      FileProperty(this,             "MARKER_POSITION_FILENAME", "", "./markerPositions.txt", false);
	public      FileProperty                   MARKERSET_FILENAME = new      FileProperty(this,                   "MARKERSET_FILENAME", "", "./data/markers.bim", false);
	public      FileProperty                MARKERLOOKUP_FILENAME = new      FileProperty(this,                "MARKERLOOKUP_FILENAME", "", "./data/markerLookup.bml", false);
	public      FileProperty                  SAMPLELIST_FILENAME = new      FileProperty(this,                  "SAMPLELIST_FILENAME", "", "./data/samples.bis", false);
	public      FileProperty               SAMPLE_SUBSET_FILENAME = new      FileProperty(this,               "SAMPLE_SUBSET_FILENAME", "", "./sampleSubset.txt", false);
	public      FileProperty                 SAMPLE_DATA_FILENAME = new      FileProperty(this,                 "SAMPLE_DATA_FILENAME", "", "./data/SampleData.txt", false);
	public      FileProperty          ORIGINAL_CENTROIDS_FILENAME = new      FileProperty(this,          "ORIGINAL_CENTROIDS_FILENAME", "", "./data/original.cent", false);
	public      FileProperty          GENOTYPE_CENTROIDS_FILENAME = new      FileProperty(this,          "GENOTYPE_CENTROIDS_FILENAME", "", "./data/genotype.cent", false);
	public      FileProperty           CHIMERA_CENTROIDS_FILENAME = new      FileProperty(this,           "CHIMERA_CENTROIDS_FILENAME", "", "./data/chimera.cent", false);
	public      FileProperty            CUSTOM_CENTROIDS_FILENAME = new      FileProperty(this,            "CUSTOM_CENTROIDS_FILENAME", "", "./data/custom.cent", false);
	public      FileProperty             DISPLAY_MARKERS_FILENAME = new      FileProperty(this,             "DISPLAY_MARKERS_FILENAME", "", "./data/test.txt", false);
	public      FileProperty            FILTERED_MARKERS_FILENAME = new      FileProperty(this,            "FILTERED_MARKERS_FILENAME", "", "./data/drops.dat", false);
	public      FileProperty                    PEDIGREE_FILENAME = new      FileProperty(this,                    "PEDIGREE_FILENAME", "", "./pedigree.dat", false);
	public      FileProperty          MOSAIC_COLOR_CODES_FILENAME = new      FileProperty(this,          "MOSAIC_COLOR_CODES_FILENAME", "", "./data/mosaic_colors.txt", false);
	public      FileProperty              MOSAIC_RESULTS_FILENAME = new      FileProperty(this,              "MOSAIC_RESULTS_FILENAME", "", "./results/Mosaicism.xln", false);
	public      FileProperty              TARGET_MARKERS_FILENAME = new      FileProperty(this,              "TARGET_MARKERS_FILENAME", "", "./targetMarkers.txt", false);
	public      FileProperty                 MOSAIC_ARMS_FILENAME = new      FileProperty(this,                 "MOSAIC_ARMS_FILENAME", "", "./MosaicArms.txt", false);
	public      FileProperty   CLUSTER_FILTER_COLLECTION_FILENAME = new      FileProperty(this,   "CLUSTER_FILTER_COLLECTION_FILENAME", "", "./data/clusterFilters.ser", false);
	public      FileProperty            SEXCHECK_RESULTS_FILENAME = new      FileProperty(this,            "SEXCHECK_RESULTS_FILENAME", "", "./results/sexCheck.xln", false);
	public      FileProperty                   GENETRACK_FILENAME = new      FileProperty(this,                   "GENETRACK_FILENAME", "", "./RefSeq.gtrack", false);
	public      FileProperty                   AB_LOOKUP_FILENAME = new      FileProperty(this,                   "AB_LOOKUP_FILENAME", "", "./AB_lookup.dat", false);
	public      FileProperty              MARKER_METRICS_FILENAME = new      FileProperty(this,              "MARKER_METRICS_FILENAME", "", "./results/markerQualityChecks.xln", false);
	public      FileProperty      MARKER_REVIEW_CRITERIA_FILENAME = new      FileProperty(this,      "MARKER_REVIEW_CRITERIA_FILENAME", "", "./results/review.criteria", false);
	public      FileProperty   MARKER_EXCLUSION_CRITERIA_FILENAME = new      FileProperty(this,   "MARKER_EXCLUSION_CRITERIA_FILENAME", "", "./results/exclusion.criteria", false);
	public      FileProperty    MARKER_COMBINED_CRITERIA_FILENAME = new      FileProperty(this,    "MARKER_COMBINED_CRITERIA_FILENAME", "", "./results/combined.criteria", false);
	public      FileProperty                  ANNOTATION_FILENAME = new      FileProperty(this,                  "ANNOTATION_FILENAME", "", "./data/annotationCollection.ser", false);
	public      FileProperty         CUSTOM_COLOR_SCHEME_FILENAME = new      FileProperty(this,         "CUSTOM_COLOR_SCHEME_FILENAME", "", "", false);
	public      FileProperty                    GC_MODEL_FILENAME = new      FileProperty(this,                    "GC_MODEL_FILENAME", "", "./data/custom.gcmodel", false);
	public      FileProperty                  COMMON_CNP_FILENAME = new      FileProperty(this,                  "COMMON_CNP_FILENAME", "", "./data/HG19 CNV edit for AGW.txt", false);
	public      FileProperty                REPORTED_CNP_FILENAME = new      FileProperty(this,                "REPORTED_CNP_FILENAME", "", "./data/HG19 Reported 2012.05.22.txt", false);
	public      FileProperty              UNREPORTED_CNP_FILENAME = new      FileProperty(this,              "UNREPORTED_CNP_FILENAME", "", "./data/HG19 Unreported 2012.05.22-2.txt", false);
	public      FileProperty                INTENSITY_PC_FILENAME = new      FileProperty(this,                "INTENSITY_PC_FILENAME", "", "./GENVISIS.PCs.extrapolated.txt", false);
	public      FileProperty                   SAMPLE_QC_FILENAME = new      FileProperty(this,                   "SAMPLE_QC_FILENAME", "", "./lrr_sd.xln", false);
	public      FileProperty          SEX_CENTROIDS_MALE_FILENAME = new      FileProperty(this,          "SEX_CENTROIDS_MALE_FILENAME", "", "", false);
	public      FileProperty        SEX_CENTROIDS_FEMALE_FILENAME = new      FileProperty(this,        "SEX_CENTROIDS_FEMALE_FILENAME", "", "", false);
	public MultiFileProperty                TWOD_LOADED_FILENAMES = new MultiFileProperty(this,                "TWOD_LOADED_FILENAMES", "", "", false);
	public MultiFileProperty                FOREST_PLOT_FILENAMES = new MultiFileProperty(this,                "FOREST_PLOT_FILENAMES", "", "", false);
	public MultiFileProperty        INDIVIDUAL_CNV_LIST_FILENAMES = new MultiFileProperty(this,        "INDIVIDUAL_CNV_LIST_FILENAMES", "", "./data/list.txt", false);
	public MultiFileProperty                REGION_LIST_FILENAMES = new MultiFileProperty(this,                "REGION_LIST_FILENAMES", "", "./data/regions.txt", false);
	public MultiFileProperty                        CNV_FILENAMES = new MultiFileProperty(this,                        "CNV_FILENAMES", "", "", false);
	public MultiFileProperty     STRATIFICATION_RESULTS_FILENAMES = new MultiFileProperty(this,     "STRATIFICATION_RESULTS_FILENAMES", "", "", false);
	public MultiFileProperty                         QQ_FILENAMES = new MultiFileProperty(this,                         "QQ_FILENAMES", "", "", false);
	public      EnumProperty                SOURCE_FILE_DELIMITER = new      EnumProperty(this,                "SOURCE_FILE_DELIMITER", "", 0, "COMMA", "TAB", "SPACE");	
	
//	private boolean jar;
	private String projectPropertiesFilename;
	private SampleList sampleList;
	private SampleData sampleData;
	private HashSet<String> cnvFilesLoadedInSampleData;
	private MarkerLookup markerLookup;
	private Logger log;
	private boolean gui;

	public Project() {
//		Files.loadProperties(this, DEFAULT_PROPERTIES, true, true, false);
		sampleList = null;
		sampleData = null;
		cnvFilesLoadedInSampleData = new HashSet<String>();
		markerLookup = null;
		log = new Logger();
		gui = false;
		this.projectPropertiesFilename = "example.properties";
	}
	
	public Project(String filename, boolean jar) {
		this(filename, null, jar);
	}
	
	// Set LOG_LEVEL to a negative value, if you do not want a log file to be generated in addition to standard out/err
	public Project(String filename, String logfile, boolean jar) {
		this();
		
		if (filename == null) {
			filename = cnv.Launch.getDefaultDebugProjectFile(true);
		}
		
		this.projectPropertiesFilename = filename;
		screenProperties();
//		Files.loadProperties(this, filename, jar, true, false);
		loadProperties(filename, jar);

//        setProperty(PROJECT_DIRECTORY, ext.verifyDirFormat(getProperty(PROJECT_DIRECTORY)));
//        setProperty(SOURCE_DIRECTORY, ext.verifyDirFormat(getProperty(SOURCE_DIRECTORY)));
//        setProperty(PROJECT_PROPERTIES_FILENAME, filename);
//		  setProperty(SAMPLE_DIRECTORY, ext.verifyDirFormat(getProperty(SAMPLE_DIRECTORY)));
        
        this.JAR_STATUS.setValue(jar);
        
		int logLevel;
		
//		logLevel = getInt(LOG_LEVEL);
		logLevel = LOG_LEVEL.getValue();
		if (logfile == null) {
			logfile = "Genvisis_"+new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log";
			if (!JAR_STATUS.getValue()) {
				logfile = PROJECT_DIRECTORY.getValue()+"logs/"+logfile;
				if (!Files.exists(PROJECT_DIRECTORY.getValue()+"logs/", JAR_STATUS.getValue())) {
					new File(PROJECT_DIRECTORY.getValue()+"logs/").mkdirs();
				}
			}
			log = new Logger(logLevel<0?null:logfile, false, Math.abs(logLevel));
		} else {
			log = new Logger(logfile, false, Math.abs(logLevel));
		}
		
	    log.report("Genvisis, v"+cnv.Launch.VERSION+"\n(c)2009-2015 Nathan Pankratz, GNU General Public License, v2\n\n"+(new Date()));
		log.report("\nCurrent project: " + getProperty(PROJECT_NAME) + "\n");
		log.report("Log level (verbosity) is set to " + getProperty(LOG_LEVEL) + "\n");
	}

	public Logger getLog() {
		return log;
	}
	
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public String getPropertyFilename() {
		return projectPropertiesFilename;
	}

	public void setPropertyFilename(String projectPropertiesFilename) {
		this.projectPropertiesFilename = projectPropertiesFilename;
	}

	public MarkerSet getMarkerSet() {
		if (Files.exists(MARKERSET_FILENAME.getValue(), JAR_STATUS.getValue())) {
			return MarkerSet.load(MARKERSET_FILENAME.getValue(), JAR_STATUS.getValue());
		} else {
			return null;
		}
	}

	public String[] getMarkerNames() {
		return getMarkerSet().getMarkerNames();
	}

	public MarkerLookup getMarkerLookup() {
		if (markerLookup == null) {
			if (Files.exists(MARKERLOOKUP_FILENAME.getValue(), JAR_STATUS.getValue())) {
				markerLookup = MarkerLookup.load(MARKERLOOKUP_FILENAME.getValue(), JAR_STATUS.getValue());
			} else {
				System.out.println("Failed to find MarkerLookup; generating one...");
				TransposeData.recreateMarkerLookup(this);
				if (Files.exists(MARKERLOOKUP_FILENAME.getValue(), JAR_STATUS.getValue())) {
					markerLookup = MarkerLookup.load(MARKERLOOKUP_FILENAME.getValue(), JAR_STATUS.getValue());
				} else {
					log.reportError("Also failed to create MarkerLookup; failing");
				}
			}
		}
		return markerLookup;
	}

	public SampleList getSampleList() {
		if (sampleList == null) {
			if (Files.exists(SAMPLELIST_FILENAME.getValue(false, false), JAR_STATUS.getValue())) {
				sampleList = SampleList.load(SAMPLELIST_FILENAME.getValue(), JAR_STATUS.getValue());
			}
			if (sampleList == null) {
				log.report("Failed to find SampleList; generating one...");
				sampleList = SampleList.generateSampleList(this);
			}
			if (sampleList != null && sampleList.getSamples().length == 0) {
				log.report("SampleList is of length zero; generating a new one...");
				sampleList = SampleList.generateSampleList(this);
			}
		}
		
		return sampleList;
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
	
	public boolean[] getSamplesToExclude() {
		boolean[] samplesToExclude;
		String[] samples;
		SampleData sampleData;
		int counter = 0;
		
		sampleData = getSampleData(0, false);
		samples = getSamples();
		samplesToExclude = new boolean[samples.length];
		for (int i = 0; i < samples.length; i++) {
			samplesToExclude[i] = sampleData.individualShouldBeExcluded(samples[i]);
			if (samplesToExclude[i]) {
				counter++;
			}
		}

		log.report("Number of samples excluded is " + counter + " (out of "+samplesToExclude.length+" total samples)");

		return samplesToExclude;
	}

	public boolean[] getSamplesToInclude(String fileWithListOfSamplesToUse) {
		return getSamplesToInclude(fileWithListOfSamplesToUse, true);
	}
	
	/**
	 * @param fileWithListOfSamplesToUse
	 *            set filename to null to only include samples not marked in the "Excluded" column of SampleData.txt
	 * @param verbose
	 *            report number to be included
	 */
	public boolean[] getSamplesToInclude(String fileWithListOfSamplesToUse, boolean verbose) {
		boolean[] samplesToInclude;
		String[] samples;
		SampleData sampleData;
		int counter = 0;
		HashSet<String> hash;
		
		if (fileWithListOfSamplesToUse != null) {
			hash = HashVec.loadFileToHashSet(fileWithListOfSamplesToUse, false);
		} else {
			hash = null;
		}
		
		sampleData = getSampleData(0, false);
		samples = getSamples();
		samplesToInclude = new boolean[samples.length];
		for (int i = 0; i < samples.length; i++) {
			if (hash == null) {
				samplesToInclude[i] = !sampleData.individualShouldBeExcluded(samples[i]);
			} else {
				samplesToInclude[i] = hash.contains(samples[i]);
			}
			if (samplesToInclude[i]) {
				counter++;
			}
		}
		
		if (verbose) {
			log.report("Number of samples to be included is " + counter + " (out of " + samplesToInclude.length + " total samples)");
		}
		
		return samplesToInclude;
	}
	
	public Sample getFullSampleFromRandomAccessFile(String sample) {
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, JAR_STATUS.getValue());
		} else {
			return null;
		}
	}

	public Sample getFullSampleFromSerialized(String sample) {
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true)+sample+".fsamp", JAR_STATUS.getValue())) {
			return Sample.loadFromSerialized(SAMPLE_DIRECTORY.getValue(false, true)+sample+".fsamp", JAR_STATUS.getValue());
		} else {
			return null;
		}
	}
	
	public Sample getPartialSampleFromRandomAccessFile(String sample) {
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, false, false, true, true, false, JAR_STATUS.getValue());
		} else {
			return null;
		}
	}
	
	public Sample getPartialSampleFromRandomAccessFile(String sample, boolean gc, boolean xy, boolean baf, boolean lrr, boolean geno) {
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, gc, xy, baf, lrr, geno, JAR_STATUS.getValue());
		} else {
			return null;
		}
	}
	
	public void resetSampleData() {
		sampleData = null;
	}

	public SampleData getSampleData(int numberOfBasicClassesToLoad, boolean loadCNVs) {
		return getSampleData(numberOfBasicClassesToLoad, loadCNVs ? this.CNV_FILENAMES.getValue() : null);
	}

	public SampleData getSampleData(int numberOfBasicClassesToLoad, String[] cnvFilenames) {
		if (cnvFilenames != null) {
			for (int i = 0; i < cnvFilenames.length; i++) {
				if (!cnvFilesLoadedInSampleData.contains(cnvFilenames[i])) {
					resetSampleData();
				}
			}
		}
		
		if (sampleData == null) {
			sampleData = new SampleData(this, numberOfBasicClassesToLoad, cnvFilenames);
//			System.err.println("SampleData loaded with "+(cnvFilenames == null?"no cnv files":Array.toStr(cnvFilenames, "/")));
			cnvFilesLoadedInSampleData = HashVec.loadToHashSet(cnvFilenames);
		}
		return sampleData;
	}
			
	public Hashtable<String,String> getFilteredHash() {
		if (getProperty(FILTERED_MARKERS_FILENAME).equals("")) {
			return new Hashtable<String,String>();
		} else if (Files.exists(FILTERED_MARKERS_FILENAME.getValue(), JAR_STATUS.getValue())) {
			return HashVec.loadFileToHashString(FILTERED_MARKERS_FILENAME.getValue(), 0, new int[] {0}, "", false, JAR_STATUS.getValue());
		} else {
			System.err.println("Error - '"+FILTERED_MARKERS_FILENAME.getValue(false, false)+"' not found");
			return new Hashtable<String,String>();
		}
	}


	public Vector<String> getStratResults() {
		String[] files;
		Vector<String> v;

		files = Files.list(PROJECT_DIRECTORY.getValue(), ".mds", JAR_STATUS.getValue());
	
		v = new Vector<String>();
		if (files == null) {
			System.err.println("Error - no .mds files found in directory");
		} else {
			for (int i = 0; i<files.length; i++) {
				v.add(PROJECT_DIRECTORY.getValue()+files[i]);
				System.out.println(PROJECT_DIRECTORY.getValue()+files[i]);
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
	
	public ClusterFilterCollection getClusterFilterCollection() {
		String filename;
		
		filename = this.CLUSTER_FILTER_COLLECTION_FILENAME.getValue(false, false);
        if (Files.exists(filename)) {
        	return ClusterFilterCollection.load(filename, this.JAR_STATUS.getValue());
        } else {
        	log.reportError("Warning - could not find cluster filter file, assuming no markers have been reclustered ("+filename+")");
        	return null;
        }
	}
	
	public AnnotationCollection getAnnotationCollection() {
		String filename;
		
		filename = this.ANNOTATION_FILENAME.getValue();
        if (Files.exists(filename)) {
    		System.out.println("Loading annotation from: "+filename);
        	return AnnotationCollection.load(filename, this.JAR_STATUS.getValue());
        } else {
        	return null;
        }
	}
	
	public String getNameOfProject() {
		return getProperty(this.PROJECT_NAME);
	}

	public String[] getIndividualRegionLists() {
		return INDIVIDUAL_CNV_LIST_FILENAMES.getValue(false);
	}
	
	public String[] getRegionLists() {
		return REGION_LIST_FILENAMES.getValue(false);
	}
	
	public String[] getPropertyKeys() {
		ArrayList<String> propList = new ArrayList<String>();
		for (Field f : this.getClass().getFields()) {
			if (containsKey(f.getName())) {
				propList.add(f.getName());
			}
		}
		return propList.toArray(new String[propList.size()]);
	}
	
	public void screenProperties() {
		BufferedReader reader;
        String trav, key;
        boolean changed;
        Vector<String> preknowns, unknowns, corrections;
        int index;
        
        changed = false;
//        knowns = Array.toStringVector(HashVec.getKeys(this, false, false));
        preknowns = Array.toStringVector(getPropertyKeys());
        unknowns = new Vector<String>();
        corrections = new Vector<String>();
        try {
	        reader = new BufferedReader(new FileReader(projectPropertiesFilename));
	        while (reader.ready()) {
	        	trav = reader.readLine();
	        	index = trav.trim().indexOf("=");
	        	if (trav.startsWith("#") || trav.startsWith("??_") || trav.equals("")) {
	        		corrections.add(trav);
	        		if (index > 0 && containsKey(trav.substring(1, index))) {
	        			preknowns.remove(trav.substring(1, index));
	        		}
	        	} else if (index > 0) {
	        		key = trav.trim().substring(0, index);

	        		if (getProperty(key) == null) {
	        			unknowns.add(key);
	        			corrections.add("??_"+trav);
	        		} else {
	        			preknowns.remove(key);
	        			corrections.add(trav);
	        		}
	        	} else {
    				log.reportError("Error - invalid property in "+projectPropertiesFilename+": "+trav);
    				log.reportError("        there is no equals sign, and comments need to be preceeded by a #");
    				log.reportError("        creating a new file with this line commented out");
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+projectPropertiesFilename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+projectPropertiesFilename+"\"");
	        System.exit(2);
        }
        
        if (unknowns.size() > 0) {
        	log.reportError("Error - check spelling for the following unexpected propert"+(unknowns.size()==1?"y":"ies")+" in "+projectPropertiesFilename+":");
        	changed = true;
        }
        for (int i = 0; i < unknowns.size(); i++) {
        	log.reportError("        "+unknowns.elementAt(i));
		}
        
        if (preknowns.size() > 0) {
        	changed = true;
        	corrections.add("");
        	corrections.add("# A few more parameters that were not originally defined:");
            for (int i = 0; i < preknowns.size(); i++) {
            	corrections.add("#"+preknowns.elementAt(i)+"="+getProperty(preknowns.elementAt(i)).getDefaultValueString());
    		}        
        }
        
        if (changed) {
            Files.backup(ext.removeDirectoryInfo(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename)+"backup/", true);
        	Files.writeList(Array.toStringArray(corrections), projectPropertiesFilename);
        }
	}

	public void saveProperties() {
		saveProperties(getPropertyFilename());
	}
	
	public void saveProperties(String outfile) {
		BufferedReader reader;
        String trav;
        boolean changed;
        Vector<String> loaded, props, changes;
        String key;
        int index;
        
        props = new Vector<String>();
        changes = new Vector<String>();
//        loaded = Array.toStringVector(HashVec.getKeys(this, false, false));
        loaded = Array.toStringVector(getPropertyKeys());
        loaded.remove(PROJECT_PROPERTIES_FILENAME);
        try {
	        reader = new BufferedReader(new FileReader(projectPropertiesFilename));
	        while (reader.ready()) {
	        	trav = reader.readLine();
	        	index = trav.trim().indexOf("=");
	        	if (trav.startsWith("#") || trav.startsWith("??_") || trav.equals("")) {
	        		props.add(trav);
	        	} else if (index > 0) {
	        		key = trav.trim().substring(0, index);
//	        		if (getProperty(key) == null) {
        			if (!containsKey(key)) {
	        			log.reportError("Unknown property '"+trav+"' not caught at startup");
        			} else {
        				props.add(key+"="+getProperty(key).getValueString());
	        			loaded.remove(key);
	        			if (!getProperty(key).getValueString().equals(trav.trim().substring(index+1))) {
	        				changes.add(key+"="+getProperty(key).getValueString());
	    					System.out.println("Was '"+trav.trim().substring(index+1)+"' now '"+getProperty(key).getValueString()+"'");
	        			}
	        		}
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+projectPropertiesFilename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+projectPropertiesFilename+"\"");
	        System.exit(2);
        }
        
        changed = false;
        if (loaded.size() > 0) {
//        	System.out.println("Default properties were unchanged - no properties file will be written");
//        	defaultProps = new Properties();
//    		Files.loadProperties(defaultProps, DEFAULT_PROPERTIES, true, true, false);
    		for (int i = 0; i < loaded.size(); i++) {
    			key = loaded.elementAt(i);
    			if (!getProperty(key).getDefaultValueString().equals(getProperty(key).getValueString())) {
    				if (!changed) {
    	            	props.add("");
    	            	props.add("# Properties where the values now differ from the defaults:");
    					changed = true;
    				}
    				props.add(key+"="+getProperty(key).getValueString());
    				changes.add(key+"="+getProperty(key).getValueString());
    				System.out.println("Default is '"+getProperty(key).getDefaultValueString()+"' now '"+getProperty(key).getValueString()+"'");
    			}
			}
        }

        if (changes.size() > 0) {
        	log.report("Changes were made to the following propert"+(changes.size()==1?"y":"ies")+" in "+projectPropertiesFilename+":");
            for (int i = 0; i < changes.size(); i++) {
            	log.report("        "+changes.elementAt(i));
    		}

            Files.backup(ext.removeDirectoryInfo(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename)+"backup/", outfile.equals(projectPropertiesFilename));
        	Files.writeList(Array.toStringArray(props), outfile);
        }
	}
	

	public void loadProperties(String filename, boolean jar) {
		InputStream is;
		
		try {
			if (jar) {
//				if (verbose) {
//					System.out.println("Loading '"+filename+"'");
//				}
				is = ClassLoader.getSystemResourceAsStream(filename);
			} else {
				is = new FileInputStream(filename);
			}
			BufferedReader reader = new BufferedReader(new InputStreamReader(is));
			String line = "";
			while ((line = reader.readLine()) != null) {
				if ("".equals(line) || (line.startsWith("#") && !line.startsWith("##"))) {
					continue;
				}
				String[] parts = line.split("=");
				if (containsKey(parts[0])) {
					if (parts.length > 1) {
						setProperty(parts[0], parts[1]);
					}
				} else {
					System.err.println("Error - unknown property key found: {" + line + "}");
				}
			}
			reader.close();
        } catch (FileNotFoundException fnfe) {
//        	if (verbose) {
        		System.err.println("Error: \""+filename+"\" could not be found");
//        	}
//			if (kill) {
//				System.exit(1);
//			}
        } catch (IOException ioe) {
			System.err.println("Error - failed to load "+filename);
			ioe.printStackTrace();
//			if (kill) {
//				System.exit(1);
//			}
        }
	}
	
	
//	public void saveProperties(String outfile) {
//		BufferedReader reader;
//        String trav;
//        boolean changed;
//        Vector<String> loaded, props, changes;
//        Properties defaultProps;
//        String key;
//        int index;
//        
//        props = new Vector<String>();
//        changes = new Vector<String>();
////        loaded = Array.toStringVector(HashVec.getKeys(this, false, false));
//        loaded = Array.toStringVector(getPropertyKeys());
//        loaded.remove(PROJECT_PROPERTIES_FILENAME);
//        try {
//	        reader = new BufferedReader(new FileReader(projectPropertiesFilename));
//	        while (reader.ready()) {
//	        	trav = reader.readLine();
//	        	index = trav.trim().indexOf("=");
//	        	if (trav.startsWith("#") || trav.startsWith("??_") || trav.equals("")) {
//	        		props.add(trav);
//	        	} else if (index > 0) {
//	        		key = trav.trim().substring(0, index);
//	        		if (getProperty(key) == null) {
//	        			log.reportError("Unknown property '"+trav+"' not caught at startup");
//        			} else {
//        				props.add(key+"="+getProperty(key));
//	        			loaded.remove(key);
//	        			if (!getProperty(key).equals(trav.trim().substring(index+1))) {
//	        				changes.add(key+"="+getProperty(key));
//	    					System.out.println("Was '"+trav.trim().substring(index+1)+"' now '"+getProperty(key)+"'");
//	        			}
//	        		}
//	        	}
//	        }
//	        reader.close();
//        } catch (FileNotFoundException fnfe) {
//	        System.err.println("Error: file \""+projectPropertiesFilename+"\" not found in current directory");
//	        System.exit(1);
//        } catch (IOException ioe) {
//	        System.err.println("Error reading file \""+projectPropertiesFilename+"\"");
//	        System.exit(2);
//        }
//        
//        changed = false;
//        if (loaded.size() > 0) {
//        	defaultProps = new Properties();
//    		Files.loadProperties(defaultProps, DEFAULT_PROPERTIES, true, true, false);
//    		for (int i = 0; i < loaded.size(); i++) {
//    			key = loaded.elementAt(i);
//    			if (!getProperty(key).equals(defaultProps.getProperty(key)) && defaultProps.getProperty(key) != null) {
//    				if (!changed) {
//    	            	props.add("");
//    	            	props.add("# Properties where the values now differ from the defaults:");
//    					changed = true;
//    				}
//    				props.add(key+"="+getProperty(key));
//       				changes.add(key+"="+getProperty(key));
//					System.out.println("Default is '"+defaultProps.getProperty(key)+"' now '"+getProperty(key)+"'");
//    			}
//			}
//        }
//
//        if (changes.size() > 0) {
//        	log.report("Changes were made to the following propert"+(changes.size()==1?"y":"ies")+" in "+projectPropertiesFilename+":");
//            for (int i = 0; i < changes.size(); i++) {
//            	log.report("        "+changes.elementAt(i));
//    		}
//
//            Files.backup(ext.removeDirectoryInfo(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename), ext.parseDirectoryOfFile(projectPropertiesFilename)+"backup/", outfile.equals(projectPropertiesFilename));
//        	Files.writeList(Array.toStringArray(props), outfile);
//        }
//	}

	public String[] getTargetMarkers() {
		String targetMarkers;
		String[] targets;

		targetMarkers = this.TARGET_MARKERS_FILENAME.getValue(false, false);
		if (new File(targetMarkers).exists()) {
			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, true);
		} else {
			if (! targetMarkers.equals("")) {
				log.report("FYI, since target markers file '" + targetMarkers + "' was not found, all markers will be exported/analyzed");
			}
			targets = null;
		}

		return targets;
	}

	public String archiveFile(String filename) {
		String backup;
		
		backup = this.BACKUP_DIRECTORY.getValue(true, true)+ext.removeDirectoryInfo(filename) + "." + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()));
		new File(filename).renameTo(new File(backup));
		
		if (Files.exists(backup)) {
			return backup;
		}
		
		log.reportError("Error - failed to backup '"+filename+"' to "+backup);
		return null;
	}

	public void message(String str) {
		message(str, "Error", JOptionPane.ERROR_MESSAGE);
	}
	
	public void setGuiState(boolean state) {
		gui = state;
	}
	
	public void message(String str, String windowTitle, int messageIcon) {
		log.reportError(str);
		if (gui) {
			JOptionPane.showMessageDialog(null, str, windowTitle, messageIcon);
		}
	}
	
	public String getLocationOfSNP_Map() {
		String filename;
		
		if (Files.exists(PROJECT_DIRECTORY.getValue()+"SNP_Map.csv")) {
			filename = PROJECT_DIRECTORY.getValue()+"SNP_Map.csv";
		} else if (Files.exists(PROJECT_DIRECTORY.getValue()+"SNP_Map.csv.gz")) {
			filename = PROJECT_DIRECTORY.getValue()+"SNP_Map.csv.gz";
		} else if (Files.exists(this.SOURCE_DIRECTORY.getValue(false, false)+"SNP_Map.csv")) {
			filename = this.SOURCE_DIRECTORY.getValue(false, false)+"SNP_Map.csv";
		} else if (Files.exists(this.SOURCE_DIRECTORY.getValue(false, false)+"SNP_Map.csv.gz")) {
			filename = this.SOURCE_DIRECTORY.getValue(false, false)+"SNP_Map.csv.gz";
		} else {
			log.reportError("Failed; could not find \"SNP_Map.csv\" or \"SNP_Map.csv.gz\" in "+PROJECT_DIRECTORY.getValue()+" or in "+this.SOURCE_DIRECTORY.getValue(false, false));
			return null;
		}
		
		return filename;
	}

	/**
	 * Grab the {@link PrincipalComponentsResiduals} from {@link Project#INTENSITY_PC_FILENAME}, will return null if can not be found
	 */
	public PrincipalComponentsResiduals loadPcResids() {
//		String pcFile = getFilename(this.INTENSITY_PC_FILENAME);
		String pcFile = this.INTENSITY_PC_FILENAME.getValue();
		PrincipalComponentsResiduals pcResids;
		if (Files.exists(pcFile)) {
			//getLog().reportTimeInfo("loading Intensity PC File " + ext.removeDirectoryInfo(pcFile));
//			pcResids = new PrincipalComponentsResiduals(this, pcFile, null, Integer.parseInt(getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false, false, null);
			pcResids = new PrincipalComponentsResiduals(this, pcFile, null, this.INTENSITY_PC_NUM_COMPONENTS.getValue(), false, 0, false, false, null);
		} else {
			getLog().reportError("Warning - did not find Intensity PC File " + pcFile + " as defined by" + this.INTENSITY_PC_FILENAME.getName() + "=" + this.INTENSITY_PC_FILENAME.getValue());
			pcResids = null;
		}
		return pcResids;
	}
	
	/**
	 * Attempts to return the gene track file from the properties, and then attempts the statgen default (only valid for local)
	 * 
	 * @param report
	 *            log which file is used
	 * 
	 */
	public String getGeneTrackFileName(boolean report) {
		String tmp = this.GENETRACK_FILENAME.getValue(false, false);
		if (tmp == null || !Files.exists(tmp)) {
			tmp = GeneSet.DIRECTORY + GeneSet.REFSEQ_TRACK;
			if (tmp == null || !Files.exists(tmp)) {
				tmp = null;
			}
		}
		if (report) {
			if (tmp == null) {
				log.reportTimeWarning("Did not find a gene track to use");
			} else {
				log.reportTimeInfo("Using gene track " + tmp);
			}
		}
		return tmp;
	}
	
	/**
	 * @return Hashtable with the indices of each marker in the project
	 */
	public Hashtable<String, Integer> getMarkerIndices() {
		String[] markerNames = getMarkerNames();
		Hashtable<String, Integer> indices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			indices.put(markerNames[i], i);
		}
		return indices;
	}
}
