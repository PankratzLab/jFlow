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

public class Project /* extends Properties*/ {
//	private static final long serialVersionUID = 1L;

//	public static final String DEFAULT_PROPERTIES = "cnv/filesys/default.properties";

//	public static final String PROJECT_NAME = "PROJECT_NAME";
//	public static final String PROJECT_PROPERTIES_FILENAME = "FILENAME";
//	public static final String PROJECT_DIRECTORY = "PROJECT_DIRECTORY";
//	public static final String SOURCE_DIRECTORY = "SOURCE_DIRECTORY";
//	public static final String SOURCE_FILENAME_EXTENSION = "SOURCE_FILENAME_EXTENSION";
//	public static final String ID_HEADER = "ID_HEADER";
//	public static final String PARSE_AT_AT_SYMBOL = "PARSE_AT_AT_SYMBOL";
//	public static final String JAR_STATUS = "JAR_STATUS";
//	public static final String SAMPLE_DIRECTORY = "SAMPLE_DIRECTORY";
//	public static final String DATA_DIRECTORY = "DATA_DIRECTORY";
//	public static final String MARKER_DATA_DIRECTORY = "MARKER_DATA_DIRECTORY";
//	public static final String RESULTS_DIRECTORY = "RESULTS_DIRECTORY";
//	public static final String DEMO_DIRECTORY = "DEMO_DIRECTORY";
//	public static final String MARKER_POSITION_FILENAME = "MARKER_POSITION_FILENAME";
//	public static final String MARKERSET_FILENAME = "MARKERSET_FILENAME";
//	public static final String MARKERLOOKUP_FILENAME = "MARKERLOOKUP_FILENAME";
//	public static final String SAMPLELIST_FILENAME = "SAMPLELIST_FILENAME";
//	public static final String SAMPLE_SUBSET_FILENAME = "SAMPLE_SUBSET_FILENAME";
//	public static final String SAMPLE_DATA_FILENAME = "SAMPLE_DATA_FILENAME";
//	public static final String ORIGINAL_CENTROIDS_FILENAME = "ORIGINAL_CENTROIDS_FILENAME";
//	public static final String GENOTYPE_CENTROIDS_FILENAME = "GENOTYPE_CENTROIDS_FILENAME";
//	public static final String CHIMERA_CENTROIDS_FILENAME = "CHIMERA_CENTROIDS_FILENAME";
//	public static final String CUSTOM_CENTROIDS_FILENAME = "CUSTOM_CENTROIDS_FILENAME";
//	public static final String SEX_CENTROIDS_FILENAMES = "SEX_CENTROIDS_FILENAMES";
//	public static final String DISPLAY_MARKERS_FILENAME = "DISPLAY_MARKERS_FILENAME";
//	public static final String FILTERED_MARKERS_FILENAME = "FILTERED_MARKERS_FILENAME";
//	public static final String PEDIGREE_FILENAME = "PEDIGREE_FILENAME";
//	public static final String MOSAIC_COLOR_CODES_FILENAME = "MOSAIC_COLOR_CODES_FILENAME";
//	public static final String MOSAIC_RESULTS_FILENAME = "MOSAIC_RESULTS_FILENAME";
//	public static final String INDIVIDUAL_CNV_LIST_FILENAMES = "INDIVIDUAL_CNV_LIST_FILENAMES";
//	public static final String REGION_LIST_FILENAMES = "REGION_LIST_FILENAMES";
//	public static final String CNV_FILENAMES = "CNV_FILENAMES";
//	public static final String STRATIFICATION_RESULTS_FILENAMES = "STRATIFICATION_RESULTS_FILENAMES";
//	public static final String GC_THRESHOLD = "GC_THRESHOLD";
//	public static final String QQ_FILENAMES = "QQ_FILENAMES";
//	public static final String DISPLAY_QUANTILES = "DISPLAY_QUANTILES";
//	public static final String DISPLAY_STANDARD_QQ = "DISPLAY_STANDARD_QQ";
//	public static final String DISPLAY_ROTATED_QQ = "DISPLAY_ROTATED_QQ";
//	public static final String TARGET_MARKERS_FILENAME = "TARGET_MARKERS_FILENAME";
//	public static final String SOURCE_FILE_DELIMITER = "SOURCE_FILE_DELIMITER";
//	public static final String PENNCNV_EXECUTABLE_DIRECTORY = "PENNCNV_EXECUTABLE_DIRECTORY";
//	public static final String PENNCNV_DATA_DIRECTORY = "PENNCNV_DATA_DIRECTORY";
//	public static final String PENNCNV_GZIP_YESNO = "PENNCNV_GZIP_YESNO";
//	public static final String PENNCNV_RESULTS_DIRECTORY = "PENNCNV_RESULTS_DIRECTORY";
//	public static final String LRRSD_CUTOFF = "LRRSD_CUTOFF";
//	public static final String MOSAIC_ARMS_FILENAME = "MOSAIC_ARMS_FILENAME";
//	public static final String NUM_THREADS = "NUM_THREADS";
//	public static final String LONG_FORMAT = "LONG_FORMAT";
//	public static final String CLUSTER_FILTER_COLLECTION_FILENAME = "CLUSTER_FILTER_COLLECTION_FILENAME";
//	public static final String SEXCHECK_RESULTS_FILENAME = "SEXCHECK_RESULTS_FILENAME";
//	public static final String TWOD_LOADED_FILENAMES = "TWOD_LOADED_FILENAMES";
//	public static final String TWOD_LOADED_VARIABLES = "TWOD_LOADED_VARIABLES";
//	public static final String GENETRACK_FILENAME = "GENETRACK_FILENAME";
//	public static final String AB_LOOKUP_FILENAME = "AB_LOOKUP_FILENAME";
//	public static final String WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER = "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER";
//	public static final String MAX_MARKERS_LOADED_PER_CYCLE = "MAX_MARKERS_LOADED_PER_CYCLE";
//	public static final String MAX_MEMORY_USED_TO_LOAD_MARKER_DATA = "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA";
//	public static final String MARKER_METRICS_FILENAME = "MARKER_METRICS_FILENAME";
//	public static final String MARKER_REVIEW_CRITERIA_FILENAME = "MARKER_REVIEW_CRITERIA_FILENAME";
//	public static final String MARKER_EXCLUSION_CRITERIA_FILENAME = "MARKER_EXCLUSION_CRITERIA_FILENAME";
//	public static final String MARKER_COMBINED_CRITERIA_FILENAME = "MARKER_COMBINED_CRITERIA_FILENAME";
//	public static final String ANNOTATION_FILENAME = "ANNOTATION_FILENAME";
//	public static final String CUSTOM_COLOR_SCHEME_FILENAME = "ANNOTATION_FILENAME";
//	public static final String BACKUP_DIRECTORY = "BACKUP_DIRECTORY";
//	public static final String SHIFT_SEX_CHR_COLORS_YESNO = "SHIFT_SEX_CHR_COLORS_YESNO";
//	public static final String QQ_MAX_NEG_LOG10_PVALUE = "QQ_MAX_NEG_LOG10_PVALUE";
//	public static final String GC_MODEL_FILENAME = "GC_MODEL_FILENAME";
//	public static final String COMMON_CNP_FILENAME = "COMMON_CNP_FILENAME";
//	public static final String REPORTED_CNP_FILENAME = "REPORTED_CNP_FILENAME";
//	public static final String UNREPORTED_CNP_FILENAME = "UNREPORTED_CNP_FILENAME";
//	public static final String LOG_LEVEL = "LOG_LEVEL";
//	public static final String INTENSITY_PC_FILENAME = "INTENSITY_PC_FILENAME";
//	public static final String INTENSITY_PC_NUM_COMPONENTS = "INTENSITY_PC_NUM_COMPONENTS";
//	public static final String FOREST_PLOT_FILENAMES = "FOREST_PLOT_FILENAMES";
//	public static final String SAMPLE_QC_FILENAME = "SAMPLE_QC_FILENAME";
	
	private static String[] verifyFiles(String[] strs, boolean isDir) {
		String[] verified = new String[strs.length];
		for (int i = 0; i < strs.length; i++) {
			verified[i] = isDir ? ext.verifyDirFormat(strs[i]) : ext.replaceAllWith(strs[i], "\\", "/");
		}
		return verified;
	}
	
	public abstract static class Property<T> {
		private final String name;
		private final T defaultValue;
		private T value;
		public Property(String name, T defVal) {
			this.name = name;
			this.defaultValue = defVal;
			this.value = defaultValue;
		}
		public      T        getValue() { return value; }
		public   void setValue(T value) { this.value = value; }
		public String         getName() { return name; }
		public      T getDefaultValue() { return defaultValue; }
		public abstract void parseValue(String valueStr);
		public String getValueString() { return value.toString(); }
		public String getDefaultValueString() { return defaultValue.toString(); }
	}
	
	public static class StringProperty extends Property<String> {
		public StringProperty(String name, String defVal) {
			super(name, defVal);
		}
		@Override
		public void parseValue(String valueStr) { setValue(valueStr); }
	}
	
	public static class StringListProperty extends Property<String[]> {
		static String delim = ";";
		public StringListProperty(String name, String[] defVal) {
			super(name, defVal);
		}
		public StringListProperty(String name, String defVal) {
			super(name, defVal.split(delim));
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
		public BooleanProperty(String name, Boolean defVal) {
			super(name, defVal);
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
		public IntegerProperty(String name, int min, int max, int defValue) {
			super(name, defValue);
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
		public DoubleProperty(String name, double min, double max, double defValue) {
			super(name, defValue);
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
//		public FileProperty(String name, File defVal, boolean dirOnly) {
//			super(name, defVal);
//			isDir = dirOnly;
//		}
		public FileProperty(String name, String defVal, boolean dirOnly) {
			super(name, dirOnly ? ext.verifyDirFormat(defVal) : ext.replaceAllWith(defVal, "\\", "/")/* == null || "".equals(defVal) ? null : new File(defVal)*/);
			isDir = dirOnly;
		}
		@Override
		public void setValue(String value) {
			super.setValue(isDir ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/"));
		}
	}
	
	public static class MultiFileProperty extends Property<String[]> {
		final boolean isDir;
//		public MultiFileProperty(String name, File[] defVal, boolean dirOnly) {
//			super(name, defVal);
//			isDir = dirOnly;
//		}
		public MultiFileProperty(String name, String[] defVal, boolean dirOnly) {
			super(name, verifyFiles(defVal, dirOnly)/*createFiles(defVal)*/);
			isDir = dirOnly;
		}
		public MultiFileProperty(String name, String defVal, boolean dirOnly) {
			super(name, verifyFiles(defVal.split(";"), dirOnly) /*createFiles(defVal)*/);
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
	}
	
	public static class EnumProperty extends StringProperty {
		String[] options;
		public EnumProperty(String name, int defaultIndex, String... opts) {
			this(name, opts[defaultIndex], opts);
		}
		public EnumProperty(String name, String defVal, String... opts) {
			super(name, defVal);
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
			throw new RuntimeException(e);
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
	
	public   IntegerProperty                            LOG_LEVEL = new   IntegerProperty(                           "LOG_LEVEL", -1, 12, 1);
	public    StringProperty                         PROJECT_NAME = new    StringProperty(                        "PROJECT_NAME", "New Project");
	public    StringProperty            SOURCE_FILENAME_EXTENSION = new    StringProperty(           "SOURCE_FILENAME_EXTENSION", ".csv");
	public    StringProperty                            ID_HEADER = new    StringProperty(                           "ID_HEADER", "Sample Name");
	public    StringProperty                TWOD_LOADED_VARIABLES = new    StringProperty(               "TWOD_LOADED_VARIABLES", "");
	public    StringProperty                            FID_ALIAS = new    StringProperty(                           "FID_ALIAS", "FID;F_ID;FamID;Fam_ID;Family;FamilyID;Family_ID");
	public    StringProperty                            IID_ALIAS = new    StringProperty(                           "IID_ALIAS", "ID;IID;I_ID;IndID;Ind_ID");
	public    StringProperty                         SAMPLE_ALIAS = new    StringProperty(                        "SAMPLE_ALIAS", "Sample;DNA;DNA#");
	public   BooleanProperty                   PARSE_AT_AT_SYMBOL = new   BooleanProperty(                  "PARSE_AT_AT_SYMBOL", Boolean.FALSE);
	public   BooleanProperty                           JAR_STATUS = new   BooleanProperty(                          "JAR_STATUS", Boolean.FALSE);
	public   BooleanProperty                    DISPLAY_QUANTILES = new   BooleanProperty(                   "DISPLAY_QUANTILES", Boolean.FALSE);
	public   BooleanProperty                  DISPLAY_STANDARD_QQ = new   BooleanProperty(                 "DISPLAY_STANDARD_QQ", Boolean.TRUE);
	public   BooleanProperty                   DISPLAY_ROTATED_QQ = new   BooleanProperty(                  "DISPLAY_ROTATED_QQ", Boolean.FALSE);
	public   BooleanProperty                   PENNCNV_GZIP_YESNO = new   BooleanProperty(                  "PENNCNV_GZIP_YESNO", Boolean.TRUE);
	public   BooleanProperty                          LONG_FORMAT = new   BooleanProperty(                         "LONG_FORMAT", Boolean.FALSE);
	public   BooleanProperty           SHIFT_SEX_CHR_COLORS_YESNO = new   BooleanProperty(          "SHIFT_SEX_CHR_COLORS_YESNO", Boolean.TRUE);
	public    DoubleProperty                         GC_THRESHOLD = new    DoubleProperty(                        "GC_THRESHOLD", 0.0, 1.0, 0.15);
	public    DoubleProperty                         LRRSD_CUTOFF = new    DoubleProperty(                        "LRRSD_CUTOFF", 0.0, 1.0, 0.32);
	public   IntegerProperty                          NUM_THREADS = new   IntegerProperty(                         "NUM_THREADS", 1, 99, 1);
	public   IntegerProperty              QQ_MAX_NEG_LOG10_PVALUE = new   IntegerProperty(             "QQ_MAX_NEG_LOG10_PVALUE", 1, 10000, 100);
	public   IntegerProperty WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER = new   IntegerProperty("WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER", 1, 1000000, 10000);
	public   IntegerProperty         MAX_MARKERS_LOADED_PER_CYCLE = new   IntegerProperty(        "MAX_MARKERS_LOADED_PER_CYCLE", 1, 10000, 100);
	public   IntegerProperty  MAX_MEMORY_USED_TO_LOAD_MARKER_DATA = new   IntegerProperty( "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA", 8, 65536, 250);
	public   IntegerProperty          INTENSITY_PC_NUM_COMPONENTS = new   IntegerProperty(         "INTENSITY_PC_NUM_COMPONENTS", 0, 10000, 100);
	public      FileProperty                    PROJECT_DIRECTORY = new      FileProperty(                   "PROJECT_DIRECTORY", "./", true);
	public      FileProperty                     SOURCE_DIRECTORY = new      FileProperty(                    "SOURCE_DIRECTORY", "./", true);
	public      FileProperty                     SAMPLE_DIRECTORY = new      FileProperty(                    "SAMPLE_DIRECTORY", "samples/", true);
	public      FileProperty                       DATA_DIRECTORY = new      FileProperty(                      "DATA_DIRECTORY", "data/", true);
	public      FileProperty                MARKER_DATA_DIRECTORY = new      FileProperty(               "MARKER_DATA_DIRECTORY", "transposed/", true);
	public      FileProperty                    RESULTS_DIRECTORY = new      FileProperty(                   "RESULTS_DIRECTORY", "results/", true);
	public      FileProperty                       DEMO_DIRECTORY = new      FileProperty(                      "DEMO_DIRECTORY", "demo/", true);
	public      FileProperty         PENNCNV_EXECUTABLE_DIRECTORY = new      FileProperty(        "PENNCNV_EXECUTABLE_DIRECTORY", "/home/npankrat/bin/", true);
	public      FileProperty               PENNCNV_DATA_DIRECTORY = new      FileProperty(              "PENNCNV_DATA_DIRECTORY", "penn_data/", true);
	public      FileProperty            PENNCNV_RESULTS_DIRECTORY = new      FileProperty(           "PENNCNV_RESULTS_DIRECTORY", "penncnv/", true);
	public      FileProperty                     BACKUP_DIRECTORY = new      FileProperty(                    "BACKUP_DIRECTORY", "backup/", true);
	public      FileProperty          PROJECT_PROPERTIES_FILENAME = new      FileProperty(                            "FILENAME", "example.properties", false);
	public      FileProperty             MARKER_POSITION_FILENAME = new      FileProperty(            "MARKER_POSITION_FILENAME", "markerPositions.txt", false);
	public      FileProperty                   MARKERSET_FILENAME = new      FileProperty(                  "MARKERSET_FILENAME", "data/markers.bim", false);
	public      FileProperty                MARKERLOOKUP_FILENAME = new      FileProperty(               "MARKERLOOKUP_FILENAME", "data/markerLookup.bml", false);
	public      FileProperty                  SAMPLELIST_FILENAME = new      FileProperty(                 "SAMPLELIST_FILENAME", "data/samples.bis", false);
	public      FileProperty               SAMPLE_SUBSET_FILENAME = new      FileProperty(              "SAMPLE_SUBSET_FILENAME", "sampleSubset.txt", false);
	public      FileProperty                 SAMPLE_DATA_FILENAME = new      FileProperty(                "SAMPLE_DATA_FILENAME", "data/SampleData.txt", false);
	public      FileProperty          ORIGINAL_CENTROIDS_FILENAME = new      FileProperty(         "ORIGINAL_CENTROIDS_FILENAME", "data/original.cent", false);
	public      FileProperty          GENOTYPE_CENTROIDS_FILENAME = new      FileProperty(         "GENOTYPE_CENTROIDS_FILENAME", "data/genotype.cent", false);
	public      FileProperty           CHIMERA_CENTROIDS_FILENAME = new      FileProperty(          "CHIMERA_CENTROIDS_FILENAME", "data/chimera.cent", false);
	public      FileProperty            CUSTOM_CENTROIDS_FILENAME = new      FileProperty(           "CUSTOM_CENTROIDS_FILENAME", "data/custom.cent", false);
	public      FileProperty             DISPLAY_MARKERS_FILENAME = new      FileProperty(            "DISPLAY_MARKERS_FILENAME", "data/test.txt", false);
	public      FileProperty            FILTERED_MARKERS_FILENAME = new      FileProperty(           "FILTERED_MARKERS_FILENAME", "data/drops.dat", false);
	public      FileProperty                    PEDIGREE_FILENAME = new      FileProperty(                   "PEDIGREE_FILENAME", "pedigree.dat", false);
	public      FileProperty          MOSAIC_COLOR_CODES_FILENAME = new      FileProperty(         "MOSAIC_COLOR_CODES_FILENAME", "data/mosaic_colors.txt", false);
	public      FileProperty              MOSAIC_RESULTS_FILENAME = new      FileProperty(             "MOSAIC_RESULTS_FILENAME", "results/Mosaicism.xln", false);
	public      FileProperty              TARGET_MARKERS_FILENAME = new      FileProperty(             "TARGET_MARKERS_FILENAME", "targetMarkers.txt", false);
	public      FileProperty                 MOSAIC_ARMS_FILENAME = new      FileProperty(                "MOSAIC_ARMS_FILENAME", "MosaicArms.txt", false);
	public      FileProperty   CLUSTER_FILTER_COLLECTION_FILENAME = new      FileProperty(  "CLUSTER_FILTER_COLLECTION_FILENAME", "data/clusterFilters.ser", false);
	public      FileProperty            SEXCHECK_RESULTS_FILENAME = new      FileProperty(           "SEXCHECK_RESULTS_FILENAME", "results/sexCheck.xln", false);
	public      FileProperty                   GENETRACK_FILENAME = new      FileProperty(                  "GENETRACK_FILENAME", "RefSeq.gtrack", false);
	public      FileProperty                   AB_LOOKUP_FILENAME = new      FileProperty(                  "AB_LOOKUP_FILENAME", "AB_lookup.dat", false);
	public      FileProperty              MARKER_METRICS_FILENAME = new      FileProperty(             "MARKER_METRICS_FILENAME", "results/markerQualityChecks.xln", false);
	public      FileProperty      MARKER_REVIEW_CRITERIA_FILENAME = new      FileProperty(     "MARKER_REVIEW_CRITERIA_FILENAME", "results/review.criteria", false);
	public      FileProperty   MARKER_EXCLUSION_CRITERIA_FILENAME = new      FileProperty(  "MARKER_EXCLUSION_CRITERIA_FILENAME", "results/exclusion.criteria", false);
	public      FileProperty    MARKER_COMBINED_CRITERIA_FILENAME = new      FileProperty(   "MARKER_COMBINED_CRITERIA_FILENAME", "results/combined.criteria", false);
	public      FileProperty                  ANNOTATION_FILENAME = new      FileProperty(                 "ANNOTATION_FILENAME", "data/annotationCollection.ser", false);
	public      FileProperty         CUSTOM_COLOR_SCHEME_FILENAME = new      FileProperty(        "CUSTOM_COLOR_SCHEME_FILENAME", "", false);
	public      FileProperty                    GC_MODEL_FILENAME = new      FileProperty(                   "GC_MODEL_FILENAME", "data/custom.gcmodel", false);
	public      FileProperty                  COMMON_CNP_FILENAME = new      FileProperty(                 "COMMON_CNP_FILENAME", "data/HG19 CNV edit for AGW.txt", false);
	public      FileProperty                REPORTED_CNP_FILENAME = new      FileProperty(               "REPORTED_CNP_FILENAME", "data/HG19 Reported 2012.05.22.txt", false);
	public      FileProperty              UNREPORTED_CNP_FILENAME = new      FileProperty(             "UNREPORTED_CNP_FILENAME", "data/HG19 Unreported 2012.05.22-2.txt", false);
	public      FileProperty                INTENSITY_PC_FILENAME = new      FileProperty(               "INTENSITY_PC_FILENAME", "GENVISIS.PCs.extrapolated.txt", false);
	public      FileProperty                   SAMPLE_QC_FILENAME = new      FileProperty(                  "SAMPLE_QC_FILENAME", "lrr_sd.xln", false);
	public MultiFileProperty              SEX_CENTROIDS_FILENAMES = new MultiFileProperty(             "SEX_CENTROIDS_FILENAMES", "", false);
	public MultiFileProperty                TWOD_LOADED_FILENAMES = new MultiFileProperty(               "TWOD_LOADED_FILENAMES", "", false);
	public MultiFileProperty                FOREST_PLOT_FILENAMES = new MultiFileProperty(               "FOREST_PLOT_FILENAMES", "", false);
	public MultiFileProperty        INDIVIDUAL_CNV_LIST_FILENAMES = new MultiFileProperty(       "INDIVIDUAL_CNV_LIST_FILENAMES", "data/list.txt", false);
	public MultiFileProperty                REGION_LIST_FILENAMES = new MultiFileProperty(               "REGION_LIST_FILENAMES", "data/regions.txt", false);
	public MultiFileProperty                        CNV_FILENAMES = new MultiFileProperty(                       "CNV_FILENAMES", "", false);
	public MultiFileProperty     STRATIFICATION_RESULTS_FILENAMES = new MultiFileProperty(    "STRATIFICATION_RESULTS_FILENAMES", "", false);
	public MultiFileProperty                         QQ_FILENAMES = new MultiFileProperty(                        "QQ_FILENAMES", "", false);
	public      EnumProperty                SOURCE_FILE_DELIMITER = new      EnumProperty(               "SOURCE_FILE_DELIMITER", 0, "COMMA", "TAB", "SPACE");	
	
	private boolean jar;
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
        
        this.jar = jar;
        
		int logLevel;
		
		logLevel = getInt(LOG_LEVEL);
		if (logfile == null) {
			logfile = "Genvisis_"+new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log";
			if (!getJarStatus()) {
				logfile = getProjectDir()+"logs/"+logfile;
				if (!Files.exists(getProjectDir()+"logs/", getJarStatus())) {
					new File(getProjectDir()+"logs/").mkdirs();
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
	
	public String getDir(FileProperty prop) {
		return getDir(prop, false, true);
	}
	public String getDir(FileProperty prop, boolean mkdirs) {
		return getDir(prop, mkdirs, true);
	}
	public String getDir(FileProperty prop, boolean mkdirs, boolean verbose) {
		return getDir(prop.getName(), mkdirs, verbose);
	}
	
	public String getDir(String property) {
		return getDir(property, false, true);
	}
	
	public String getDir(String property, boolean mkdirs) {
		return getDir(property, mkdirs, true);
	}
	
	public String getDir(String property, boolean mkdirs, boolean verbose) {
		String dir = null;
		
		if (containsKey(property)) {
//			dir = getProperty(property);
			dir = this.<FileProperty>getProperty(property).getValue();
			dir = ext.replaceTilde(dir);
			if (!dir.startsWith("~") && !dir.startsWith("/") && dir.indexOf(":") == -1) {
				dir = getProperty(PROJECT_DIRECTORY)+dir;
			}
			if (!Files.exists(dir, jar)) {
				if (mkdirs && !jar) {
					new File(dir).mkdirs();
				} else if (verbose) {
					log.reportError("Error - directory '"+dir+"' does not exist");
				}
			}
		} else {
			log.reportError("Error - directory property '"+property+"' is undefined in cnv.filesys.Project");
		}
		
		return dir;
	}
	
	public String getPropertyFilename() {
		return projectPropertiesFilename;
	}

	public void setPropertyFilename(String projectPropertiesFilename) {
		this.projectPropertiesFilename = projectPropertiesFilename;
	}

	public String getProjectDir() {
		return getProperty(PROJECT_DIRECTORY);
	}
	
	public String getFilename(FileProperty prop) {
		return getFilename(prop, false, true);
	}
	public String getFilename(FileProperty prop, boolean make, boolean verbose) {
		return getFilename(prop, null, make, verbose);
	}
	public String getFilename(FileProperty prop, String subdir, boolean make, boolean verbose) {
		return getFilename(prop.getName(), subdir, make, verbose);
	}

	public String getFilename(String property) {
		return getFilename(property, false, true);
	}
	
	public String getFilename(String property, boolean make, boolean verbose) {
		return getFilename(property, null, make, verbose);
	}
	
	public String getFilename(String property, String subdirectory, boolean make, boolean verbose) {
		String file = null;
		
		if (containsKey(property)) {
//			file = getProperty(property);
			file = this.<FileProperty>getProperty(property).getValue();
			if (!file.startsWith("/") && file.indexOf(":") == -1) {
				file = getProperty(PROJECT_DIRECTORY)+(subdirectory==null?"":getProperty(subdirectory))+file;
			}
			if (!Files.exists(file, getJarStatus())) {
				if (make) {
					new File(ext.parseDirectoryOfFile(file)).mkdirs();
				} else if (verbose) {
					if (!Files.exists(ext.parseDirectoryOfFile(file), getJarStatus())) {
						System.err.println("Error - the directory ('"+ext.parseDirectoryOfFile(file)+"') of the file you're trying to access/create ('"+ext.removeDirectoryInfo(file)+"') does not exist");
					} else {
						System.err.println("Error - file '"+file+"' does not exist");
					}
				}
			}
		} else {
			System.err.println("Error - file property '"+property+"' is undefined in cnv.filesys.Project");
		}

		return file;
	}
	
	public String[] getFilenames(MultiFileProperty prop) {
		return getFilenames(prop, false);
	}
	public String[] getFilenames(MultiFileProperty prop, boolean suppressMissing) {
		return getFilenames(prop.getName(), suppressMissing);
	}
	
	public String[] getFilenames(String type) {
		return getFilenames(type, false);
	}
	
	public String[] getFilenames(String type, boolean suppressMissing) {
		String[] files = null;
		Vector<String> v;
		
		v = new Vector<String>();
		if (containsKey(type)) {
//			files = getProperty(type).split(";");
			files = this.<MultiFileProperty>getProperty(type).getValue();
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
	
//	public double getDouble(String variable) {
//		String trav;
//		
//		trav = getProperty(variable);
//		try {
//			return Double.parseDouble(trav);
//		} catch (NumberFormatException nfe) {
//			System.err.println("Error - '"+trav+"' is not a valid value for "+variable);
//			return Double.NaN;
//		}		
//	}
	public double getDouble(String variable) {
		return this.<Property<? extends Number>>getProperty(variable).getValue().doubleValue();
//		try {
//			return Double.parseDouble(trav);
//		} catch (NumberFormatException nfe) {
//			System.err.println("Error - '"+trav+"' is not a valid value for "+variable);
//			return Double.NaN;
//		}		
	}
	public double getDouble(Property<? extends Number> prop) {
		return prop.getValue().doubleValue();
	}
	
//	public float getFloat(String variable) {
//		String trav;
//		
//		trav = getProperty(variable);
//		try {
//			return Float.parseFloat(trav);
//		} catch (NumberFormatException nfe) {
//			System.err.println("Error - '"+trav+"' is not a valid value for "+variable);
//			return Float.NaN;
//		}		
//	}
	public float getFloat(String variable) {
		return this.<Property<? extends Number>>getProperty(variable).getValue().floatValue();
	}
	public float getFloat(Property<? extends Number> prop) {
		return prop.getValue().floatValue();
	}

//	public int getInt(String variable) {
//		String trav;
//		
//		trav = getProperty(variable);
//		try {
//			return Integer.parseInt(trav);
//		} catch (NumberFormatException nfe) {
//			System.err.println("Error - '"+trav+"' is not a valid value for "+variable);
//			return Integer.MIN_VALUE;
//		}		
//	}
	public int getInt(String variable) {
		return this.<Property<? extends Number>>getProperty(variable).getValue().intValue();
	}
	public int getInt(Property<? extends Number> prop) {
		return prop.getValue().intValue();
	}

//	public boolean getBoolean(String variable) {
//		String trav;
//		trav = getProperty(variable);
//		if (trav.toLowerCase().equals("true") || trav.toLowerCase().equals("t") || trav.toLowerCase().equals("yes") || trav.toLowerCase().equals("y") || trav.equals("1")) {
//			return true;
//		}
//		if (trav.toLowerCase().equals("false") || trav.toLowerCase().equals("f") || trav.toLowerCase().equals("no") || trav.toLowerCase().equals("n") || trav.equals("0")) {
//			return false;
//		}
//		System.err.println("Error - '"+trav+"' is not a valid boolean value for "+variable);
//		return false;
//	}
	public boolean getBoolean(Property<Boolean> prop) {
		return prop.getValue();
	}
	public boolean getBoolean(String variable) {
		Boolean trav;
		trav = this.<Property<Boolean>>getProperty(variable).getValue();
		return trav;
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
		if (markerLookup == null) {
			if (Files.exists(getFilename(MARKERLOOKUP_FILENAME), getJarStatus())) {
				markerLookup = MarkerLookup.load(getFilename(MARKERLOOKUP_FILENAME), getJarStatus());
			} else {
				System.out.println("Failed to find MarkerLookup; generating one...");
				TransposeData.recreateMarkerLookup(this);
				if (Files.exists(getFilename(MARKERLOOKUP_FILENAME), getJarStatus())) {
					markerLookup = MarkerLookup.load(getFilename(MARKERLOOKUP_FILENAME), getJarStatus());
				} else {
					log.reportError("Also failed to create MarkerLookup; failing");
				}
			}
		}
		return markerLookup;
	}

	public SampleList getSampleList() {
		if (sampleList == null) {
			if (Files.exists(getFilename(SAMPLELIST_FILENAME, false, false), getJarStatus())) {
				sampleList = SampleList.load(getFilename(SAMPLELIST_FILENAME), getJarStatus());
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
		if (Files.exists(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, getJarStatus())) {
			return Sample.loadFromRandomAccessFile(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, getJarStatus());
		} else {
			return null;
		}
	}

	public Sample getFullSampleFromSerialized(String sample) {
		if (Files.exists(getDir(SAMPLE_DIRECTORY)+sample+".fsamp", getJarStatus())) {
			return Sample.loadFromSerialized(getDir(SAMPLE_DIRECTORY)+sample+".fsamp", getJarStatus());
		} else {
			return null;
		}
	}
	
	public Sample getPartialSampleFromRandomAccessFile(String sample) {
		if (Files.exists(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, getJarStatus())) {
			return Sample.loadFromRandomAccessFile(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, false, false, true, true, false, getJarStatus());
		} else {
			return null;
		}
	}
	
	public Sample getPartialSampleFromRandomAccessFile(String sample, boolean gc, boolean xy, boolean baf, boolean lrr, boolean geno) {
		if (Files.exists(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, getJarStatus())) {
			return Sample.loadFromRandomAccessFile(getDir(SAMPLE_DIRECTORY) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, gc, xy, baf, lrr, geno, getJarStatus());
		} else {
			return null;
		}
	}
	
	public void resetSampleData() {
		sampleData = null;
	}

	public SampleData getSampleData(int numberOfBasicClassesToLoad, boolean loadCNVs) {
		return getSampleData(numberOfBasicClassesToLoad, loadCNVs ? getFilenames(this.CNV_FILENAMES) : null);
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
	
	public ClusterFilterCollection getClusterFilterCollection() {
		String filename;
		
		filename = getFilename(this.CLUSTER_FILTER_COLLECTION_FILENAME, false, false);
        if (Files.exists(filename)) {
        	return ClusterFilterCollection.load(filename, jar);
        } else {
        	log.reportError("Warning - could not find cluster filter file, assuming no markers have been reclustered ("+filename+")");
        	return null;
        }
	}
	
	public AnnotationCollection getAnnotationCollection() {
		String filename;
		
		filename = getFilename(this.ANNOTATION_FILENAME);
        if (Files.exists(filename)) {
    		System.out.println("Loading annotation from: "+filename);
        	return AnnotationCollection.load(filename, jar);
        } else {
        	return null;
        }
	}
	
	public String getNameOfProject() {
		return getProperty(this.PROJECT_NAME);
	}

	public boolean getJarStatus() {
		return jar;
	}
	
	public void setJarStatus(boolean jar) {
		this.jar = jar;
	}
	
	public String[] getIndividualRegionLists() {
		return getFilenames(INDIVIDUAL_CNV_LIST_FILENAMES, false);
	}
	
	public String[] getRegionLists() {
		return getFilenames(REGION_LIST_FILENAMES, false);
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
        Vector<String> knowns, unknowns, corrections;
        int index;
        
        changed = false;
//        knowns = Array.toStringVector(HashVec.getKeys(this, false, false));
        knowns = Array.toStringVector(getPropertyKeys());
        unknowns = new Vector<String>();
        corrections = new Vector<String>();
        try {
	        reader = new BufferedReader(new FileReader(projectPropertiesFilename));
	        while (reader.ready()) {
	        	trav = reader.readLine();
	        	index = trav.trim().indexOf("=");
	        	if (trav.startsWith("#") || trav.startsWith("??_") || trav.equals("")) {
	        		corrections.add(trav);
	        		if (index > 0 && getProperty(trav.substring(1, index)) != null) {
	        			knowns.remove(trav.substring(1, index));
	        		}
	        	} else if (index > 0) {
	        		key = trav.trim().substring(0, index);

	        		if (getProperty(key) == null) {
	        			unknowns.add(key);
	        			corrections.add("??_"+trav);
	        		} else {
	        			knowns.remove(key);
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
        
        if (knowns.size() > 0) {
        	changed = true;
        	corrections.add("");
        	corrections.add("# A few more parameters that were not originally defined:");
            for (int i = 0; i < knowns.size(); i++) {
            	corrections.add("#"+knowns.elementAt(i)+"="+getProperty(knowns.elementAt(i)).getDefaultValueString());
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
    				props.add(key+"="+getProperty(key).getValue());
       				changes.add(key+"="+getProperty(key).getValue());
					System.out.println("Default is '"+getProperty(key).getDefaultValue()+"' now '"+getProperty(key).getValue()+"'");
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
					setProperty(parts[0], parts[1]);
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

		targetMarkers = getFilename(this.TARGET_MARKERS_FILENAME, false, false);
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
		
		backup = getDir(this.BACKUP_DIRECTORY, true, true)+ext.removeDirectoryInfo(filename) + "." + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()));
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
		
		if (Files.exists(getProjectDir()+"SNP_Map.csv")) {
			filename = getProjectDir()+"SNP_Map.csv";
		} else if (Files.exists(getProjectDir()+"SNP_Map.csv.gz")) {
			filename = getProjectDir()+"SNP_Map.csv.gz";
		} else if (Files.exists(getDir(this.SOURCE_DIRECTORY, false, false)+"SNP_Map.csv")) {
			filename = getDir(this.SOURCE_DIRECTORY, false, false)+"SNP_Map.csv";
		} else if (Files.exists(getDir(this.SOURCE_DIRECTORY, false, false)+"SNP_Map.csv.gz")) {
			filename = getDir(this.SOURCE_DIRECTORY, false, false)+"SNP_Map.csv.gz";
		} else {
			log.reportError("Failed; could not find \"SNP_Map.csv\" or \"SNP_Map.csv.gz\" in "+getProjectDir()+" or in "+getDir(this.SOURCE_DIRECTORY, false, false));
			return null;
		}
		
		return filename;
	}

	/**
	 * Grab the {@link PrincipalComponentsResiduals} from {@link Project#INTENSITY_PC_FILENAME}, will return null if can not be found
	 */
	public PrincipalComponentsResiduals loadPcResids() {
		String pcFile = getFilename(this.INTENSITY_PC_FILENAME);
		PrincipalComponentsResiduals pcResids;
		if (Files.exists(pcFile)) {
			//getLog().reportTimeInfo("loading Intensity PC File " + ext.removeDirectoryInfo(pcFile));
//			pcResids = new PrincipalComponentsResiduals(this, pcFile, null, Integer.parseInt(getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false, false, null);
			pcResids = new PrincipalComponentsResiduals(this, pcFile, null, getProperty(this.INTENSITY_PC_NUM_COMPONENTS), false, 0, false, false, null);
		} else {
			getLog().reportError("Warning - did not find Intensity PC File " + pcFile + " as defined by" + this.INTENSITY_PC_FILENAME + "=" + this.getProperty(this.INTENSITY_PC_FILENAME));
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
		String tmp = getFilename(this.GENETRACK_FILENAME, false, false);
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
