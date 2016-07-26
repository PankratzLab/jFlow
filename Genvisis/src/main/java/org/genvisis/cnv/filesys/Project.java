package org.genvisis.cnv.filesys;

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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map.Entry;
import java.util.Vector;

import javax.swing.JOptionPane;
import javax.swing.JProgressBar;

import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.CurrentManifest;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneSet;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;

public class Project {

//	private static String[] verifyFiles(String[] strs, boolean isDir) {
//		String[] verified = new String[strs.length];
//		for (int i = 0; i < strs.length; i++) {
//			verified[i] = isDir ? ext.verifyDirFormat(strs[i]) : ext.replaceAllWith(strs[i], "\\", "/");
//		}
//		return verified;
//	}
	
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
		@Override
		public String toString() { return getName() + "=" + getValueString(); }
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
		boolean isFile = false;
		boolean isDir = false;
		public StringListProperty(Project proj, String name, String desc, String[] defVal, boolean file, boolean dir) {
			super(proj, name, desc, defVal);
			isFile = file;
			isDir = dir;
		}
		public StringListProperty(Project proj, String name, String desc, String defVal, boolean file, boolean dir) {
			super(proj, name, desc, defVal.equals("") ? new String[0] : defVal.split(delim));
			isFile = file;
			isDir = dir;
		}
		@Override
		public void parseValue(String valueStr) {
		    String[] pts = "".equals(valueStr) ? new String[0] : valueStr.split(delim);
			this.setValue(pts);
		}
		@Override
		public String getValueString() {
			return Array.toStr(getValue(), delim);
		}
		public String getDefaultValueString() { return Array.toStr(getDefaultValue(), delim); }
		public String getDelimiter() { return delim; }
		@Override
		public String[] getValue() {
		    String[] values = super.getValue();
		    if (isFile || isDir) {
    		    for (int i = 0; i < values.length; i++) {
    		        if (!"".equals(values[i]) && !values[i].startsWith(".") && !values[i].startsWith("/") && values[i].indexOf(":") == -1) {
    		            values[i] = getProject().PROJECT_DIRECTORY.getValue() + values[i];
    	            }
    		    }
		    }
		    return values;
		}
		@Override
		public void setValue(String[] value) {
		    String[] values = value;
		    if (isFile || isDir) {
                for (int i = 0; i < values.length; i++) {
                    if (values[i].startsWith(getProject().PROJECT_DIRECTORY.getValue())) {
                        values[i] = values[i].substring(getProject().PROJECT_DIRECTORY.getValue().length());
                    }
                }
            }
		    super.setValue(values);
		}
		public void removeValue(String value) {
		    String[] newValues = new String[this.getValue().length - 1];
		    String[] values = this.getValue();
		    int index = 0;
		    for (int i = 0; i < values.length; i++) {
		        boolean skip = false;
		        if ((isFile || isDir) && (ext.verifyDirFormat(values[i]).equals(ext.verifyDirFormat(value)))) { 
		            skip = true;
		        } else if (values[i].equals(value)) {
		            skip = true;
		        }
		        if (skip) {
		            continue;
		        }
		        newValues[index++] = values[i];
		    }
		    this.setValue(newValues);
		}
		public void addValue(String value) {
		    this.setValue(Array.addStrToArray(value, getValue(), 0));
		}
		public void addValue(String valu, int index) {
            String value = valu;
            if (isDir) {
                value = ext.verifyDirFormat(value);
            } else if (isFile) {
                value = ext.verifyDirFormat(value);
                value = value.substring(0, value.length() - 1);
            }
		    this.setValue(Array.addStrToArray(value, getValue(), index));
		}
	}
	
	public static class BooleanProperty extends Property<Boolean> {
		public BooleanProperty(Project proj, String name, String desc, Boolean defVal) {
			super(proj, name, desc, defVal);
		}
		@Override
		public void parseValue(String valueStr) { 
			Boolean newValue = valueStr.equals("") ? this.getDefaultValue() : Boolean.valueOf(valueStr);
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
			Integer newValue = valueStr.equals("") ? this.getDefaultValue() : Integer.valueOf(valueStr);
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
			Double newValue = valueStr.equals("") ? this.getDefaultValue() : Double.valueOf(valueStr);
			setValue(newValue);
		}
		@Override
		public void setValue(Double value) {
			if (value < min || value > max) {
				throw new RuntimeException("Error - values for property " + getName() + " must be within " + min + "-" + max + "; " + value + " is not valid");
			}
			super.setValue(value);
		}
		
		public double getMinValue(){
			return min;
		}
		
		public double getMaxValue(){
			return max;
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
		// TODO NP asks: When is this subdir option ever used?
		public String getValue(String subdir, boolean mkdirs, boolean verbose) {
			String valu = super.getValue();
			if (valu.contains("~")) {
				valu = ext.replaceTilde(valu);
			}
			
			String tempValue = valu;
			if (!"".equals(valu) && !valu.startsWith(".") && !valu.startsWith("/") && valu.indexOf(":") == -1) {
				if (isDir) {
					if (this.getName().equals("PROJECT_DIRECTORY")) { // happens with example.properties
					 	tempValue = LaunchProperties.directoryOfLaunchProperties(LaunchProperties.DEFAULT_PROPERTIES_FILE) + valu;
					} else {
						tempValue = getProject().PROJECT_DIRECTORY.getValue() + valu;
					}
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
						if (!Files.exists(ext.parseDirectoryOfFile(tempValue), getProject().JAR_STATUS.getValue())) {
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
	
	public static class EnumProperty<T extends Enum<T>> extends Property<T> {
	    T[] enumValues;
	    int defaultIndex;
	    public EnumProperty(Project proj, String name, String description, int defaultIndex, Class<T> opts) {
            super(proj, name, description, opts.getEnumConstants()[defaultIndex]);
            this.enumValues = opts.getEnumConstants();
            this.defaultIndex = defaultIndex;
        }
	    public EnumProperty(Project proj, String name, String description, T defaultOpt, Class<T> opts) {
	        super(proj, name, description, defaultOpt);
	        this.enumValues = opts.getEnumConstants();
	        for (int i = 0; i < enumValues.length; i++) {
	            if (defaultOpt == enumValues[i]) {
	                this.defaultIndex = i;
	                break;
	            }
	        }
        }

        @Override
        public void parseValue(String valueStr) {
            for (T val : enumValues) {
                if (val.toString().equals(valueStr)) {
                    this.setValue(val);
                }
            }
        }
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
	
	public   IntegerProperty                            LOG_LEVEL = new    IntegerProperty(this,                            "LOG_LEVEL", "", -1, 12, 1);
	public    StringProperty                         PROJECT_NAME = new     StringProperty(this,                         "PROJECT_NAME", "Project Name", "New Project");
	public    StringProperty            SOURCE_FILENAME_EXTENSION = new     StringProperty(this,            "SOURCE_FILENAME_EXTENSION", "", ".csv");
	public    StringProperty                            ID_HEADER = new     StringProperty(this,                            "ID_HEADER", "", "Sample Name");
	public    StringProperty                TWOD_LOADED_VARIABLES = new     StringProperty(this,                "TWOD_LOADED_VARIABLES", "", "");
	public    StringProperty                            FID_ALIAS = new     StringProperty(this,                            "FID_ALIAS", "", "FID;F_ID;FamID;Fam_ID;Family;FamilyID;Family_ID");
	public    StringProperty                            IID_ALIAS = new     StringProperty(this,                            "IID_ALIAS", "", "ID;IID;I_ID;IndID;Ind_ID");
	public    StringProperty                         SAMPLE_ALIAS = new     StringProperty(this,                         "SAMPLE_ALIAS", "", "Sample;DNA;DNA#");
	public   BooleanProperty                   PARSE_AT_AT_SYMBOL = new    BooleanProperty(this,                   "PARSE_AT_AT_SYMBOL", "", Boolean.FALSE);
	public   BooleanProperty                           JAR_STATUS = new    BooleanProperty(this,                           "JAR_STATUS", "", Boolean.FALSE);
	public   BooleanProperty                    DISPLAY_QUANTILES = new    BooleanProperty(this,                    "DISPLAY_QUANTILES", "", Boolean.FALSE);
	public   BooleanProperty                  DISPLAY_STANDARD_QQ = new    BooleanProperty(this,                  "DISPLAY_STANDARD_QQ", "", Boolean.TRUE);
	public   BooleanProperty                   DISPLAY_ROTATED_QQ = new    BooleanProperty(this,                   "DISPLAY_ROTATED_QQ", "", Boolean.FALSE);
	public   BooleanProperty                   PENNCNV_GZIP_YESNO = new    BooleanProperty(this,                   "PENNCNV_GZIP_YESNO", "", Boolean.TRUE);
	public   BooleanProperty                          LONG_FORMAT = new    BooleanProperty(this,                          "LONG_FORMAT", "", Boolean.FALSE);
	public   BooleanProperty           SHIFT_SEX_CHR_COLORS_YESNO = new    BooleanProperty(this,           "SHIFT_SEX_CHR_COLORS_YESNO", "", Boolean.TRUE);
	public    DoubleProperty        BLAST_PROPORTION_MATCH_FILTER = new     DoubleProperty(this,        "BLAST_PROPORTION_MATCH_FILTER", "", 0.0, 1.0, 0.80);
	public    DoubleProperty                         GC_THRESHOLD = new     DoubleProperty(this,                         "GC_THRESHOLD", "", 0.0, 1.0, 0.15);
	public    DoubleProperty                      XY_SCALE_FACTOR = new     DoubleProperty(this,                         "XY_SCALE_FACTOR", "", 0.001, Double.MAX_VALUE, 1);
	public    DoubleProperty                         LRRSD_CUTOFF = new     DoubleProperty(this,                         "LRRSD_CUTOFF", "", 0.0, 3.0, 0.32);
	public    DoubleProperty            SAMPLE_CALLRATE_THRESHOLD = new     DoubleProperty(this,            "SAMPLE_CALLRATE_THRESHOLD", "", 0.0, 1.0, 0.95);
	public   IntegerProperty                          NUM_THREADS = new    IntegerProperty(this,                          "NUM_THREADS", "", 1, 99, 1);
	public   IntegerProperty              QQ_MAX_NEG_LOG10_PVALUE = new    IntegerProperty(this,              "QQ_MAX_NEG_LOG10_PVALUE", "", 1, 10000, 100);
	public   IntegerProperty WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER = new    IntegerProperty(this, "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER", "", 1, 1000000, 10000);
	public   IntegerProperty         MAX_MARKERS_LOADED_PER_CYCLE = new    IntegerProperty(this,         "MAX_MARKERS_LOADED_PER_CYCLE", "", 1, 10000, 100);
	public   IntegerProperty  MAX_MEMORY_USED_TO_LOAD_MARKER_DATA = new    IntegerProperty(this,  "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA", "", 8, 65536, 250);
	public   IntegerProperty          INTENSITY_PC_NUM_COMPONENTS = new    IntegerProperty(this,          "INTENSITY_PC_NUM_COMPONENTS", "", 0, 10000, 100);
	public      FileProperty                    PROJECT_DIRECTORY = new       FileProperty(this,                    "PROJECT_DIRECTORY", "", "./", true);
	public      FileProperty                     SOURCE_DIRECTORY = new       FileProperty(this,                     "SOURCE_DIRECTORY", "", "./", true);
	public      FileProperty                     SAMPLE_DIRECTORY = new       FileProperty(this,                     "SAMPLE_DIRECTORY", "", "samples/", true);
	public      FileProperty                       DATA_DIRECTORY = new       FileProperty(this,                       "DATA_DIRECTORY", "", "data/", true);
	public      FileProperty                MARKER_DATA_DIRECTORY = new       FileProperty(this,                "MARKER_DATA_DIRECTORY", "", "transposed/", true);
	public      FileProperty                    RESULTS_DIRECTORY = new       FileProperty(this,                    "RESULTS_DIRECTORY", "", "results/", true);
	public      FileProperty                       DEMO_DIRECTORY = new       FileProperty(this,                       "DEMO_DIRECTORY", "", "demo/", true);
	public      FileProperty         PENNCNV_EXECUTABLE_DIRECTORY = new       FileProperty(this,         "PENNCNV_EXECUTABLE_DIRECTORY", "", "/home/npankrat/bin/", true);
	public      FileProperty               PENNCNV_DATA_DIRECTORY = new       FileProperty(this,               "PENNCNV_DATA_DIRECTORY", "", "penn_data/", true);
	public      FileProperty            PENNCNV_RESULTS_DIRECTORY = new       FileProperty(this,            "PENNCNV_RESULTS_DIRECTORY", "", "penncnv/", true);
	public      FileProperty                     BACKUP_DIRECTORY = new       FileProperty(this,                     "BACKUP_DIRECTORY", "", "backup/", true);
	public      FileProperty          PROJECT_PROPERTIES_FILENAME = new       FileProperty(this,                             "FILENAME", "", "example.properties", false);
	public      FileProperty             MARKER_POSITION_FILENAME = new       FileProperty(this,             "MARKER_POSITION_FILENAME", "", "markerPositions.txt", false);
	public      FileProperty                   MARKERSET_FILENAME = new       FileProperty(this,                   "MARKERSET_FILENAME", "", "data/markers.ser", false);
	public      FileProperty                MARKERLOOKUP_FILENAME = new       FileProperty(this,                "MARKERLOOKUP_FILENAME", "", "data/markerLookup.ser", false);
	public      FileProperty                  SAMPLELIST_FILENAME = new       FileProperty(this,                  "SAMPLELIST_FILENAME", "", "data/samples.ser", false);
	public      FileProperty               SAMPLE_SUBSET_FILENAME = new       FileProperty(this,               "SAMPLE_SUBSET_FILENAME", "", "sampleSubset.txt", false);
	public      FileProperty                 SAMPLE_DATA_FILENAME = new       FileProperty(this,                 "SAMPLE_DATA_FILENAME", "", "data/SampleData.txt", false);
	public      FileProperty          ORIGINAL_CENTROIDS_FILENAME = new       FileProperty(this,          "ORIGINAL_CENTROIDS_FILENAME", "", "data/original.cent", false);
	public      FileProperty          GENOTYPE_CENTROIDS_FILENAME = new       FileProperty(this,          "GENOTYPE_CENTROIDS_FILENAME", "", "data/genotype.cent", false);
	public      FileProperty           CHIMERA_CENTROIDS_FILENAME = new       FileProperty(this,           "CHIMERA_CENTROIDS_FILENAME", "", "data/chimera.cent", false);
	public      FileProperty            CUSTOM_CENTROIDS_FILENAME = new       FileProperty(this,            "CUSTOM_CENTROIDS_FILENAME", "", "data/custom.cent", false);
	public      FileProperty            FILTERED_MARKERS_FILENAME = new       FileProperty(this,            "FILTERED_MARKERS_FILENAME", "", "data/drops.dat", false);
	public      FileProperty                    PEDIGREE_FILENAME = new       FileProperty(this,                    "PEDIGREE_FILENAME", "", "pedigree.dat", false);
	public      FileProperty          MOSAIC_COLOR_CODES_FILENAME = new       FileProperty(this,          "MOSAIC_COLOR_CODES_FILENAME", "", "data/mosaic_colors.txt", false);
	public      FileProperty              MOSAIC_RESULTS_FILENAME = new       FileProperty(this,              "MOSAIC_RESULTS_FILENAME", "", "results/Mosaicism.xln", false);
	public      FileProperty   CLUSTER_FILTER_COLLECTION_FILENAME = new       FileProperty(this,   "CLUSTER_FILTER_COLLECTION_FILENAME", "", "data/clusterFilters.ser", false);
	public      FileProperty            SEXCHECK_RESULTS_FILENAME = new       FileProperty(this,            "SEXCHECK_RESULTS_FILENAME", "", "results/sexCheck.xln", false);
	public      FileProperty                   GENETRACK_FILENAME = new       FileProperty(this,                   "GENETRACK_FILENAME", "", "RefSeq.gtrack", false);
	public      FileProperty                   AB_LOOKUP_FILENAME = new       FileProperty(this,                   "AB_LOOKUP_FILENAME", "", "AB_lookup.dat", false);
	public      FileProperty              MARKER_METRICS_FILENAME = new       FileProperty(this,              "MARKER_METRICS_FILENAME", "", "results/markerQualityChecks.xln", false);
	public      FileProperty      MARKER_REVIEW_CRITERIA_FILENAME = new       FileProperty(this,      "MARKER_REVIEW_CRITERIA_FILENAME", "", "results/review.criteria", false);
	public      FileProperty   MARKER_EXCLUSION_CRITERIA_FILENAME = new       FileProperty(this,   "MARKER_EXCLUSION_CRITERIA_FILENAME", "", "results/exclusion.criteria", false);
	public      FileProperty    MARKER_COMBINED_CRITERIA_FILENAME = new       FileProperty(this,    "MARKER_COMBINED_CRITERIA_FILENAME", "", "results/combined.criteria", false);
	public      FileProperty                  ANNOTATION_FILENAME = new       FileProperty(this,                  "ANNOTATION_FILENAME", "", "data/annotationCollection.ser", false);
	public      FileProperty            BLAST_ANNOTATION_FILENAME = new       FileProperty(this,            "BLAST_ANNOTATION_FILENAME", "", "data/blast.vcf.gz", false);
	public      FileProperty         CUSTOM_COLOR_SCHEME_FILENAME = new       FileProperty(this,         "CUSTOM_COLOR_SCHEME_FILENAME", "", "", false);
	public      FileProperty                    GC_MODEL_FILENAME = new       FileProperty(this,                    "GC_MODEL_FILENAME", "", "data/custom.gcmodel", false);
	public      FileProperty                  COMMON_CNP_FILENAME = new       FileProperty(this,                  "COMMON_CNP_FILENAME", "", "data/HG19 CNV edit for AGW.txt", false);
	public      FileProperty                REPORTED_CNP_FILENAME = new       FileProperty(this,                "REPORTED_CNP_FILENAME", "", "data/HG19 Reported 2012.05.22.txt", false);
	public      FileProperty              UNREPORTED_CNP_FILENAME = new       FileProperty(this,              "UNREPORTED_CNP_FILENAME", "", "data/HG19 Unreported 2012.05.22-2.txt", false);
	public      FileProperty                INTENSITY_PC_FILENAME = new       FileProperty(this,                "INTENSITY_PC_FILENAME", "", "PCA_GENVISIS.PCs.extrapolated.txt", false);
	public      FileProperty                   SAMPLE_QC_FILENAME = new       FileProperty(this,                   "SAMPLE_QC_FILENAME", "", "lrr_sd.xln", false);
	public      FileProperty          SEX_CENTROIDS_MALE_FILENAME = new       FileProperty(this,          "SEX_CENTROIDS_MALE_FILENAME", "", "", false);
	public      FileProperty        SEX_CENTROIDS_FEMALE_FILENAME = new       FileProperty(this,        "SEX_CENTROIDS_FEMALE_FILENAME", "", "", false);
	public      FileProperty      REFERENCE_GENOME_FASTA_FILENAME = new       FileProperty(this,      "REFERENCE_GENOME_FASTA_FILENAME", "", "hg19_canonical.fa", false);
	public      FileProperty              GENOME_CLUSTER_FILENAME = new       FileProperty(this,              "GENOME_CLUSTER_FILENAME", "", "cluster.genome.gz", false); 
	public      FileProperty                  CUSTOM_PFB_FILENAME = new       FileProperty(this,                  "CUSTOM_PFB_FILENAME", "", "data/custom.pfb", false);
	public      FileProperty                         HMM_FILENAME = new       FileProperty(this,                         "HMM_FILENAME", "", "data/hhall.hmm", false);
    public      FileProperty        INTENSITY_PC_MARKERS_FILENAME = new       FileProperty(this,        "INTENSITY_PC_MARKERS_FILENAME", "", "GENVISIS.PCs.markers.txt", false);
    public StringListProperty                 GENE_LIST_FILENAMES = new StringListProperty(this,                  "GENE_LIST_FILENAMES", "", "data/genes.txt", true, false);
	public StringListProperty            TARGET_MARKERS_FILENAMES = new StringListProperty(this,             "TARGET_MARKERS_FILENAMES", "", "targetMarkers.txt", true, false);
	public StringListProperty           DISPLAY_MARKERS_FILENAMES = new StringListProperty(this,            "DISPLAY_MARKERS_FILENAMES", "", "data/test.txt", true, false);
	public StringListProperty               TWOD_LOADED_FILENAMES = new StringListProperty(this,                "TWOD_LOADED_FILENAMES", "", "", true, false);
	public StringListProperty               FOREST_PLOT_FILENAMES = new StringListProperty(this,                "FOREST_PLOT_FILENAMES", "", "", true, false);
	public StringListProperty       INDIVIDUAL_CNV_LIST_FILENAMES = new StringListProperty(this,        "INDIVIDUAL_CNV_LIST_FILENAMES", "", "data/list.txt", true, false);
	public StringListProperty               REGION_LIST_FILENAMES = new StringListProperty(this,                "REGION_LIST_FILENAMES", "", "data/regions.txt", true, false);
	public StringListProperty                       CNV_FILENAMES = new StringListProperty(this,                        "CNV_FILENAMES", "", "", true, false);
	public StringListProperty    STRATIFICATION_RESULTS_FILENAMES = new StringListProperty(this,     "STRATIFICATION_RESULTS_FILENAMES", "", "", true, false);
	public StringListProperty                        QQ_FILENAMES = new StringListProperty(this,                         "QQ_FILENAMES", "", "", true, false);
	public StringListProperty  GC_CORRECTION_PARAMETERS_FILENAMES = new StringListProperty(this,   "GC_CORRECTION_PARAMETERS_FILENAMES", "", "", true, false);
	public StringListProperty                 PLINK_DIR_FILEROOTS = new StringListProperty(this,                  "PLINK_DIR_FILEROOTS", "", "", true, false);
	public StringListProperty          MARKER_COLOR_KEY_FILENAMES = new StringListProperty(this,           "MARKER_COLOR_KEY_FILENAMES", "", "", true, false);

	public EnumProperty<SOURCE_FILE_DELIMITERS>   SOURCE_FILE_DELIMITER = new EnumProperty<SOURCE_FILE_DELIMITERS>(this, "SOURCE_FILE_DELIMITER", "", 0, SOURCE_FILE_DELIMITERS.class);	
	public EnumProperty<ARRAY>                               ARRAY_TYPE = new EnumProperty<ARRAY>(this, "ARRAY_TYPE", "", 0, ARRAY.class);	
	public EnumProperty<GENOME_BUILD>              GENOME_BUILD_VERSION = new EnumProperty<GENOME_BUILD>(this, "GENOME_BUILD_VERSION", "The build version of the genome, options are " + Arrays.asList(GENOME_BUILD.values()).toString(), 0, GENOME_BUILD.class);

	private String projectPropertiesFilename;
	private SampleList sampleList;
	private SampleData sampleData;
	private HashSet<String> cnvFilesLoadedInSampleData;
	private HashMap<String, SourceFileHeaderData> sourceFileHeaders;
	private MarkerLookup markerLookup;
	private Logger log;
	private boolean gui;
	private ProgressMonitor progressMonitor;
	public ProgressMonitor getProgressMonitor() {
	    return progressMonitor;
	}

	
	public static final String HEADERS_FILENAME = "source_headers.ser";
	
	public Project() {
		sampleList = null;
		sampleData = null;
		cnvFilesLoadedInSampleData = new HashSet<String>();
		markerLookup = null;
		log = new Logger();
		gui = false;
		this.projectPropertiesFilename = "example.properties";
		initializeProgressMonitor(null);
	}
	
	public Project(String filename, boolean jar) {
		this(filename, null, jar);
	}

	public Project(String filename, String logfile, boolean jar) {
		this(filename, logfile, jar, true);
	}
	
	// Set LOG_LEVEL to a negative value, if you do not want a log file to be generated in addition to standard out/err
	public Project(String filename, String logfile, boolean jar, boolean createHeaders) {
		this();
		
		if (filename == null) {
			filename = org.genvisis.cnv.Launch.getDefaultDebugProjectFile(true);
		}
		
		this.projectPropertiesFilename = filename;
		screenProperties();
		loadProperties(filename, jar);

//        setProperty(PROJECT_DIRECTORY, ext.verifyDirFormat(getProperty(PROJECT_DIRECTORY)));
//        setProperty(SOURCE_DIRECTORY, ext.verifyDirFormat(getProperty(SOURCE_DIRECTORY)));
//        setProperty(PROJECT_PROPERTIES_FILENAME, filename);
//		  setProperty(SAMPLE_DIRECTORY, ext.verifyDirFormat(getProperty(SAMPLE_DIRECTORY)));
        
		
        this.JAR_STATUS.setValue(jar);
        
		int logLevel;
		
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

		if (Files.exists(SAMPLE_DIRECTORY.getValue()) && (new File(SAMPLE_DIRECTORY.getValue()).list().length > 0)) {
			// skip source file headers, sample files already parsed
		} else if (createHeaders && Files.list(SOURCE_DIRECTORY.getValue(), SOURCE_FILENAME_EXTENSION.getValue(), false).length > 0) {
			HashMap<String, SourceFileHeaderData> headers = readHeadersFile(false);
			setSourceFileHeaders(headers);
		}

		
		log.report(CurrentManifest.getGenvisisInfo());
		log.report("\nJava version: " + System.getProperty("java.version"));

		try {
			if (!Files.checkJVMUpToDateApprox()) {
				log.reportError("\nYOUR VERSION OF JAVA IS OUT OF DATE; update if you get a NoSuchMethodError");
			}
		} catch (Exception e) {
			log.reportError("\nCould not parse Java version and check for possible compatibility issues\n");
		}
		
		log.report("\nCurrent project: " + getProperty(PROJECT_NAME) + "\n");
		log.report("Log level (verbosity) is set to " + getProperty(LOG_LEVEL) + "\n");
		
		
		updateProject(this);
	}

    private static void updateProject(Project proj) {
        updateProperty(proj.SAMPLELIST_FILENAME, ".bis", "sample list");
        updateProperty(proj.MARKERLOOKUP_FILENAME, ".bml", "marker lookup");
        updateProperty(proj.MARKERSET_FILENAME, ".bim", "marker set");
        proj.saveProperties(new Property[]{proj.SAMPLELIST_FILENAME, proj.MARKERLOOKUP_FILENAME, proj.MARKERSET_FILENAME});
    }
    
    private static void updateProperty(FileProperty prop, String prevExt, String fileDescriptor) {
        String file, newFile, warning = null, error = null;
        
        file = prop.getValue();
        if (!file.endsWith(".ser")) {
            newFile = ext.rootOf(file, false) + ".ser";
            if (Files.exists(file)) {
                if (!file.endsWith(prevExt)) {
                    // unlikely, but error
                    warning = "Warning - found " + fileDescriptor + " file, but with an unexpected extension.  Renaming to .ser from \"" + file + "\".";
                }
                if (Files.exists(newFile) && (new File(newFile)).length() > 0) {
                    warning = "Warning - found .ser version of " + fileDescriptor + " file; altering property value to point to .ser file at \"" + newFile + "\".";
                } else {
                    (new File(newFile)).delete();
                    (new File(file)).renameTo(new File(newFile));
                }
                prop.setValue(ext.rootOf(file, false) + ".ser");
            } else { 
                if (Files.exists(newFile)) {
                    // property set incorrectly, file exists as .ser but prop set to .bis, rename property
                    prop.setValue(ext.rootOf(file, false) + ".ser");
                } else {
                    // error, no sample list file (and property set to old value)
                    error = "Error - " + fileDescriptor + " file not found!  Set to: \"" + file + "\".";
                }
            }
        }
        if (warning != null) {
            System.err.println(warning);
        }
        if (error != null) {
            System.err.println(error);
        }
    }
	
	
	
	private boolean reasonableCheckForParsedSource() {
        String sampleDirectory = SAMPLE_DIRECTORY.getValue(false, false);
        // TODO strict check for #files == #samples?
        return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION, false).length > 0 && getSampleList() != null && getSampleList().getSamples().length > 0;
	}
	
	@SuppressWarnings("unchecked")
    private HashMap<String, SourceFileHeaderData> readHeadersFile(boolean waitIfMissing) {
	    String file = PROJECT_DIRECTORY.getValue() + "source.headers";
	    
	    if (Files.exists(file)) {
	        log.report("Found source.headers file, renaming to " + HEADERS_FILENAME);
	        (new File(file)).renameTo(new File(PROJECT_DIRECTORY.getValue() + HEADERS_FILENAME));
	    }
	    file = PROJECT_DIRECTORY.getValue() + HEADERS_FILENAME;
	    
	    if (Files.exists(file)) {
	        HashMap<String, SourceFileHeaderData> headers = (HashMap<String, SourceFileHeaderData>) SerializedFiles.readSerial(file, JAR_STATUS.getValue().booleanValue(), getLog(), false);
	        if (headers != null) {
	            return headers;
	        } else {
	            // error reading headers; let's delete
	            getLog().reportError(ext.getTime() + "]\tError reading source file header metadata.  Deleting file and reparsing.");
	            getLog().reportError(ext.getTime() + "]\tThis is only relevant if desired data columns are non-default AND source files are not yet parsed into " + Sample.SAMPLE_FILE_EXTENSION + " files.");
	            getLog().reportError(ext.getTime() + "]\tA quick check (which may be incorrect) suggest this " + (reasonableCheckForParsedSource() ? "IS LIKELY NOT " : "IS LIKELY") + " to be an issue.");
	            (new File(file)).delete();
	        }
	    } 
        if (!waitIfMissing) {
	        new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        log.report("Parsing source file headers in background thread.");
                        setSourceFileHeaders(SourceFileHeaderData.validate(SOURCE_DIRECTORY.getValue(), SOURCE_FILENAME_EXTENSION.getValue(), true, log, null));
                        log.report("Source file header parsing complete.");
                    } catch (Exception e) {
                        log.reportException(e);
                    }
                }
            }).start();
	        return null;
        } else {
            try {
                log.report("Parsing source file headers in active thread.");
                setSourceFileHeaders(SourceFileHeaderData.validate(SOURCE_DIRECTORY.getValue(), SOURCE_FILENAME_EXTENSION.getValue(), true, log, null));
                log.report("Source file header parsing complete.");
                return getSourceFileHeaders(false);
            } catch (Exception e) {
                log.reportException(e);
                return null;
            }
        }
    }
	
	private void writeHeadersFile() {
	    SerializedFiles.writeSerial(sourceFileHeaders, PROJECT_DIRECTORY.getValue() + HEADERS_FILENAME);
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
			getLog().reportFileNotFound(MARKERSET_FILENAME.getValue());
			return null;
		}
	}

	public String[] getMarkerNames() {
		return getMarkerSet().getMarkerNames();
	}

	public synchronized MarkerLookup getMarkerLookup() {
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

	public synchronized SampleList getSampleList() {
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
	
	public boolean[] getSamplesToInclude(String fileWithListOfSamplesToUse, boolean verbose) {
		return getSamplesToInclude(fileWithListOfSamplesToUse, false, verbose);
	}

	/**
	 * @param fileWithListOfSamplesToUse
	 *            set filename to null to only include samples not marked in the "Excluded" column of SampleData.txt
	 * @param overlapExclude
	 *            if a file is provided, the union of the file and samples not marked in the "Excluded" can be obtained
	 * @param verbose
	 *            report number to be included
	 */
	public boolean[] getSamplesToInclude(String fileWithListOfSamplesToUse, boolean overlapExclude, boolean verbose) {
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
			} else if (hash != null && overlapExclude) {
				samplesToInclude[i] = !sampleData.individualShouldBeExcluded(samples[i]) && hash.contains(samples[i]);
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
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, JAR_STATUS.getValue());
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
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, false, false, true, true, false, JAR_STATUS.getValue());
		} else {
			return null;
		}
	}
	
	public Sample getPartialSampleFromRandomAccessFile(String sample, boolean gc, boolean xy, boolean baf, boolean lrr, boolean geno) {
		if (Files.exists(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, JAR_STATUS.getValue())) {
			return Sample.loadFromRandomAccessFile(SAMPLE_DIRECTORY.getValue(false, true) + sample + Sample.SAMPLE_FILE_EXTENSION, gc, xy, baf, lrr, geno, JAR_STATUS.getValue());
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
					break;
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
	
	public void saveProperties(String outfile, String[] propsToSave) {
        BufferedReader reader;
        String trav;
        boolean changed;
        HashSet<String> propKeysOfInterest;
        Vector<String> loaded, props, changes;
        String key;
        int index;
        
        propKeysOfInterest = new HashSet<String>();
        for (String prop : propsToSave) {
            propKeysOfInterest.add(prop);
        }
        props = new Vector<String>();
        changes = new Vector<String>();
        loaded = Array.toStringVector(propsToSave);
        loaded.remove(PROJECT_PROPERTIES_FILENAME);
        try {
            reader = new BufferedReader(new FileReader(projectPropertiesFilename));
            trav = null;
            while ((trav = reader.readLine()) != null) {
                index = trav.trim().indexOf("=");
                if (trav.startsWith("#") || trav.startsWith("??_") || trav.equals("")) {
                    props.add(trav);
                } else if (index > 0) {
                    key = trav.trim().substring(0, index);
//                  if (getProperty(key) == null) {
                    if (!containsKey(key)) {
                        log.reportError("Unknown property '"+trav+"' not caught at startup");
                    } else if (propKeysOfInterest.contains(key)) {
                        String valueString = getProperty(key).getValueString();
                        if (!key.equals(PROJECT_DIRECTORY.getName()) && !key.equals(SOURCE_DIRECTORY.getName())) {
                            valueString = valueString.replace(PROJECT_DIRECTORY.getValue(), "");
                        }
                        props.add(key+"="+valueString);
                        loaded.remove(key);
                        if (!valueString.equals(trav.trim().substring(index+1))) {
                            changes.add(key+"="+valueString);
                            log.report("Was '"+trav.trim().substring(index+1)+"' now '"+valueString+"'");
                        }
                    } else {
                        props.add(trav);
                    }
                }
            }
            reader.close();
        } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \""+projectPropertiesFilename+"\" not found in current directory");
            System.exit(1);
        } catch (IOException ioe) {
            log.reportError("Error reading file \""+projectPropertiesFilename+"\"");
            System.exit(2);
        }
        
        changed = false;
        if (loaded.size() > 0) {
            for (int i = 0; i < loaded.size(); i++) {
                key = loaded.elementAt(i);
                
                String valueString = getProperty(key).getValueString();
                String defaultValueString = getProperty(key).getDefaultValueString();
                if (!key.equals(PROJECT_DIRECTORY.getName()) && !key.equals(SOURCE_DIRECTORY.getName())) {
                    valueString = valueString.replace(PROJECT_DIRECTORY.getValue(), "");
                    defaultValueString = defaultValueString.replace(PROJECT_DIRECTORY.getValue(), "");
                }
                
                if (!defaultValueString.equals(valueString)) {
                    if (!changed) {
                        props.add("");
                        props.add("# Properties where the values now differ from the defaults:");
                        changed = true;
                    }
                    props.add(key+"="+valueString);
                    changes.add(key+"="+valueString);
                    log.report("Default for Project property " + key + " was '"+defaultValueString+"' and is now '"+valueString+"'");
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
	
	@SuppressWarnings("rawtypes")
    public <T extends Property> void saveProperties(T[] propsToSave) {
	    String[] propNames = new String[propsToSave.length];
	    for (int i = 0; i < propsToSave.length; i++) {
	        propNames[i] = propsToSave[i].getName();
	    }
	    saveProperties(getPropertyFilename(), propNames);
	}
	
	public void saveProperties(String[] propNamesToSave) {
	    saveProperties(getPropertyFilename(), propNamesToSave);
	}

	public void saveProperties(String outfile) {
	    saveProperties(outfile, getPropertyKeys());
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
					setProperty(parts[0], parts.length > 1 ? parts[1] : "");
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
	    return getTargetMarkers(this.TARGET_MARKERS_FILENAMES.getValue()[0]);
	}
	
	public String[] getTargetMarkers(String targetMarkerFile) {
		String targetMarkers;
		String[] targets;
		
		if (targetMarkerFile == null) {
			return getMarkerNames();
		}

		targetMarkers = targetMarkerFile;
		if (new File(targetMarkers).exists()) {
			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, true);
		} else {
			if (! targetMarkers.equals("")) {
				log.report("FYI, since target markers file '" + targetMarkers + "' was not found, all markers will be exported/analyzed");
			}
			targets = getMarkerNames();
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

	public void setGuiState(boolean state) {
		gui = state;
	}
	
	public void initializeProgressMonitor(JProgressBar progBar) {
	    this.progressMonitor = new ProgressMonitor(progBar, this.log);
	}

	/**
	 *	Reports message to the log and if and only if a GUI is being used, it also creates a message dialog as well
	 * 
	 * @param str			The message to display
	 * @param windowTitle	Title of the message
	 * @param messageIcon	Icon to use can be any of the following:
	 * 									JOptionPane.ERROR_MESSAGE
	 * 									JOptionPane.INFORMATION_MESSAGE 
	 * 									JOptionPane.WARNING_MESSAGE
	 * 									JOptionPane.QUESTION_MESSAGE
	 * 									JOptionPane.PLAIN_MESSAGE
	 */
	public void message(String str, String windowTitle, int messageIcon) {
		log.reportError(str);
		if (gui) {
			JOptionPane.showMessageDialog(null, str, windowTitle, messageIcon);
		}
	}
	
	/**
	 *	Reports message to the log and if and only if a GUI is being used, it also creates a message dialog as well
	 *	This simplified method assumes this is an error and says as much in the window title it creates  
	 * 
	 * @param str			The message to display
	 * 
	*/
	public void message(String str) {
		message(str, "Error", JOptionPane.ERROR_MESSAGE);
	}

	public String getLocationOfSNP_Map(boolean verbose) {
		String filename;
		
		String projDir = PROJECT_DIRECTORY.getValue();
        String snpMap = "SNP_Map.csv";
        String snpMapGz = "SNP_Map.csv.gz";
        if (Files.exists(projDir + snpMap)) {
			filename = projDir + snpMap;
		} else if (Files.exists(projDir + snpMapGz)) {
			filename = projDir + snpMapGz;
		} else {
            String srcDir = this.SOURCE_DIRECTORY.getValue();
            if (Files.exists(srcDir + snpMap)) {
            	filename = srcDir + snpMap;
            } else if (Files.exists(srcDir + snpMapGz)) {
            	filename = srcDir + snpMapGz;
            } else {
            	if (verbose) {
            	    log.reportError("Failed; could not find \"" + snpMap + "\" or \"" + snpMapGz + "\" in " + projDir + " or in " + srcDir);
            	}
            	return null;
            }
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
			getLog().reportError("Warning - did not find Intensity PC File " + pcFile + " as defined by " + this.INTENSITY_PC_FILENAME.getName() + "=" + this.INTENSITY_PC_FILENAME.getValue());
			pcResids = null;
		}
		return pcResids;
	}

	/**
	 * @return the {@link Pedigree} if the {@link Project#PEDIGREE_FILENAME} exists, null otherwise
	 */
	public Pedigree loadPedigree() {
		String ped = PEDIGREE_FILENAME.getValue();
		Pedigree pedigree = null;
		if (!Files.exists(ped)) {
			log.reportTimeWarning("Did not find pedigree file " + ped);
		} else if ((new File(ped)).length() == 0) {
		    log.reportTimeWarning("Pedigree file " + ped + " was empty.");
	    } else {
			pedigree = new Pedigree(this); // will load from project
		}
		return pedigree;
	}
	
	/**
	 * Attempts to return the gene track file from the properties, and then attempts other default locations, set's property if found elsewhere
	 * 
	 * @param verbose	whether to report
	 * @return			GeneTrack, if it found one, otherwise null
	 */
	public String getGeneTrackFilename(boolean verbose) {
		String geneTrackFilename = this.GENETRACK_FILENAME.getValue(false, false);
		if (geneTrackFilename == null || !Files.exists(geneTrackFilename)) {
			geneTrackFilename = Files.firstPathToFileThatExists(Aliases.REFERENCE_FOLDERS, GeneSet.REFSEQ_TRACK, true, false, log);
			if (geneTrackFilename == null || !Files.exists(geneTrackFilename)) {
				geneTrackFilename = null;
			}
		}
		if (verbose) {
			if (geneTrackFilename == null) {
				log.reportTimeWarning("Did not find a GeneTrack to use");
			} else {
				log.reportTimeInfo("Using gene track " + geneTrackFilename);
			}
		}
		if (geneTrackFilename != null && Files.exists(geneTrackFilename) && !geneTrackFilename.equals(this.GENETRACK_FILENAME.getValue(false, false))) {
		    this.GENETRACK_FILENAME.setValue(geneTrackFilename);
		}
		return geneTrackFilename;
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
	
	/**
	 * @return the oulier hash from all samples
	 * @throws Exception
	 */
	public Hashtable<String, Float> loadOutliersFromSamples() throws Exception {
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		String[] samples = getSamples();
		for (int i = 0; i < samples.length; i++) {
			Hashtable<String, Float> sOutliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(SAMPLE_DIRECTORY.getValue() + samples[i] + Sample.SAMPLE_FILE_EXTENSION);
			if (sOutliers != null && sOutliers.size() > 0) {
				outliers.putAll(sOutliers);
			}
		}
		return outliers;
	}

	public String[] getNonCNMarkers() {
		String[] mkrs = getMarkerNames();
		ARRAY myArrayType = this.ARRAY_TYPE.getValue();
		ArrayList<String> nonCNs = new ArrayList<String>();
		for (int i = 0; i < mkrs.length; i++) {
			if (!myArrayType.isCNOnly(mkrs[i])) {
				nonCNs.add(mkrs[i]);
			}
		}
		return Array.toStringArray(nonCNs);
	}

	public String[] getAutosomalNonCNMarkers() {
		String[] mkrs = getAutosomalMarkers();
		ARRAY myArrayType = this.ARRAY_TYPE.getValue();
		ArrayList<String> nonCNs = new ArrayList<String>();
		for (int i = 0; i < mkrs.length; i++) {
			if (!myArrayType.isCNOnly(mkrs[i])) {
				nonCNs.add(mkrs[i]);
			}
		}
		return Array.toStringArray(nonCNs);
	}

	public boolean[] getCNMarkers() {
	    String[] mkrs = getMarkerNames();
        boolean[] cnB = new boolean[mkrs.length];
        ARRAY myArrayType = this.ARRAY_TYPE.getValue();
	    for (int i = 0; i < mkrs.length; i++) {
	        cnB[i] = myArrayType.isCNOnly(mkrs[i]);
	    }
	    return cnB;
	}
	
	public String[] getAutosomalMarkers() {
		MarkerSet markerSet = getMarkerSet();
		byte[] chrs = markerSet.getChrs();
		ArrayList<String> tmp = new ArrayList<String>();
		for (int i = 0; i < chrs.length; i++) {
			if (chrs[i] < 23 && chrs[i] > 0) {
				tmp.add(markerSet.getMarkerNames()[i]);
			}
		}
		return tmp.toArray(new String[tmp.size()]);
	}
	
	
	/**
	 * @return indices of autosomal markers
	 */
	public int[] getAutosomalMarkerIndices() {
		String[] autosomalMarkers = getAutosomalMarkers();
		int[] indices = ext.indexLargeFactors(autosomalMarkers, getMarkerNames(), true, log, true, false);
		return indices;
	}
	
	/**
	 * @return boolean representation of autosomal markers
	 */
	public boolean[] getAutosomalMarkerBoolean() {
		int[] indices = getAutosomalMarkerIndices();
		boolean[] autoB = Array.booleanArray(getMarkerNames().length, false);
		for (int i = 0; i < indices.length; i++) {
			autoB[indices[i]] = true;
		}
		return autoB;
	}
	
	
	/**
	 * For copying an existing project to a new project that will have the same essential data
	 */
	public void copyBasicFiles(Project projectToCopyTo, boolean overwrite) {
		HashSet<FileProperty> propsToCop = new HashSet<FileProperty>();
		propsToCop.add(MARKERSET_FILENAME);
		propsToCop.add(MARKER_POSITION_FILENAME);
		propsToCop.add(SAMPLE_DATA_FILENAME);
		propsToCop.add(SAMPLELIST_FILENAME);
		propsToCop.add(GC_MODEL_FILENAME);
		propsToCop.add(PEDIGREE_FILENAME);
		for (FileProperty fileProperty : propsToCop) {
			copyToNewProject(this, projectToCopyTo, fileProperty.getName(), overwrite);
		}
		// TODO copy targetMarkers files?  All?  Just default file if exists?
	}

	/**
	 * @param projOriginal
	 * @param projectToCopyTo
	 * @param fileProperty
	 *            the file property to copy
	 * @param overwrite
	 *            whether to overwrite this file if it exists in the new destination
	 * @return whether the file was copied
	 */
	private static boolean copyToNewProject(Project projOriginal, Project projectToCopyTo, String fileProperty, boolean overwrite) {
		String fileOriginal = (String) projOriginal.getProperty(fileProperty).getValue();
		String fileToCopyTo = (String) projectToCopyTo.getProperty(fileProperty).getValue();

		if (Files.exists(fileOriginal) && (!Files.exists(fileToCopyTo) || overwrite)) {
			String dir = ext.parseDirectoryOfFile(fileToCopyTo);
			new File(dir).mkdirs();
			return Files.copyFileUsingFileChannels(fileOriginal, fileToCopyTo, projOriginal.getLog());
		} else {
			return false;
		}
	}

	public enum SOURCE_FILE_DELIMITERS {
	    COMMA("[\\s]*,[\\s]*", ","),
	    TAB("[ ]*\t[ ]*", "\t"),
	    SPACE("[\\s]+", " ");
	    
	    String delim;
	    HashSet<String> alts = new HashSet<String>();
	    private SOURCE_FILE_DELIMITERS(String... delimValues) {
	        this.delim = delimValues[0];
	        for (String d : delimValues) {
	            alts.add(d);
	        }
        }
	    public String getDelimiter() {
	        return this.delim;
	    }
	    public static SOURCE_FILE_DELIMITERS getDelimiter(String value) {
	        for (SOURCE_FILE_DELIMITERS delim : SOURCE_FILE_DELIMITERS.values()) {
	            if (delim.getDelimiter().equals(value) || delim.getDelimiter() == value || delim.alts.contains(value)) {
	                return delim;
	            }
	        }
	        return null;
	    }
	}
	
	public ARRAY getArrayType() {
		return ARRAY_TYPE.getValue();
	}
	/**
	 * Used for determining how to import, compute LRR, etc
	 *
	 */
	public enum ARRAY {

		/**
		 * Your friendly Illumina arrays
		 */
		ILLUMINA(new String[] { "cnvi" }, 50/*, 650000*/), /**
		 * Supports CHP format
		 */
		AFFY_GW6(new String[] { "CN_" }, 25/*, 650000*/),
		/**
		 * Supports CHP and CNCHP formated input
		 */
		AFFY_GW6_CN(new String[] { "CN_" }, 25/*, 909622*/),

		/**
		 * For bamFiles
		 */
		NGS(new String[] { "*" }, 100/*, 0*/),
		
//		DBGAP(new String[] {}, 0, 909622)
		;

		/**
		 * Used for copy number only probe-set identification
		 */
		private String[] cnFlags;
		/**
		 * Length of the probe sequences on the array
		 */
		private int probeLength;

		private ARRAY(String[] cnFlags, int probeLength) {
			this.cnFlags = cnFlags;
			this.probeLength = probeLength;
		}

		public String[] getCnFlags() {

			return cnFlags;
		}

		public void setCnFlags(String[] cnFlags) {
			this.cnFlags = cnFlags;
		}

		public int getProbeLength() {
			return probeLength;
		}

		public void setProbeLength(int probeLength) {
			this.probeLength = probeLength;
		}

		/**
		 * @param markerName
		 * @return whether the marker trips the {@link ARRAY#cnFlags} flags
		 */
		public boolean isCNOnly(String markerName) {
			if (this == NGS) {
				return NGS_MARKER_TYPE.getType(markerName) != NGS_MARKER_TYPE.VARIANT_SITE;// only non cn type we have
			} else {
				for (int i = 0; i < cnFlags.length; i++) {
					if (ext.indexOfStartsWith(cnFlags[i], new String[] { markerName }, false) >= 0) {
						return true;
					}
				}
				return false;
			}
		}

	}

	/**
	 * Even if there are not outliers we still enforce "outliers.ser" to exist.
	 */
	public void verifyAndGenerateOutliers(boolean verbose) {
		Sample.verifyAndGenerateOutliers(this, NUM_THREADS.getValue(), false);
	}

	public HashMap<String, SourceFileHeaderData> getSourceFileHeaders(boolean readIfNull) {
	    return sourceFileHeaders == null ? readIfNull ? readHeadersFile(true) : null : sourceFileHeaders;
	}
	
    public void setSourceFileHeaders(HashMap<String, SourceFileHeaderData> sourceFileHeaders) {
        this.sourceFileHeaders = sourceFileHeaders;
        writeHeadersFile();
	}

	/**
	 * 
	 * Mainly for development of methods involving changes to sample/transposed files
	 * 
	 * @param projOriginal
	 * @param tag
	 *            this tag will serve as the new directory under the original project directory
	 * @return
	 */
	public static Project prepareNewProject(Project projOriginal, String tag) {
		String newProjectFile = ext.addToRoot(projOriginal.PROJECT_PROPERTIES_FILENAME.getValue(), "." + tag);
		Files.copyFileUsingFileChannels(projOriginal.PROJECT_PROPERTIES_FILENAME.getValue(), newProjectFile, projOriginal.getLog());
		Project projCorrected = new Project(newProjectFile, false);
		String newDir = projOriginal.PROJECT_DIRECTORY.getValue() + tag + "/";
		projOriginal.getLog().reportTimeInfo("Preparing project " + newProjectFile + " in " + newDir);
		new File(newDir).mkdirs();
		projCorrected.PROJECT_DIRECTORY.setValue(newDir);
		projCorrected.PROJECT_NAME.setValue(projOriginal.PROJECT_NAME.getValue() + tag);
		projCorrected.saveProperties();
		projOriginal.copyBasicFiles(projCorrected, false);
		return projCorrected;
	}
    
	public static void main(String[] args) {
        int numArgs = args.length;
        String filename = null;
        Project proj;
        HashMap<String, String> kvPairs = new HashMap<String, String>();
        
        String usage = "\n"+
            "cnv.filesys.Project requires 2+ arguments\n"+
            "   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
            "   (2+) key-value pairs for properties (i.e. NUM_THREADS=6 (not the default))\n"+
            "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help")
                    || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("proj=")) {
                filename = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].contains("=")) {
                String[] parts = args[i].split("=");
                if (parts.length > 2) {
                    break;
                }
                kvPairs.put(parts[0], parts.length > 1 ? parts[1] : "");
                numArgs--;
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        
        proj = new Project(filename, false);
        for (Entry<String, String> kv : kvPairs.entrySet()) {
            try {
                proj.setProperty(kv.getKey(), kv.getValue());
            } catch (Throwable e) {
                System.err.println("Error - malformed key-value property: {" + kv.getKey() + "=" + kv.getValue() + "}");
            }
        }
        proj.saveProperties();
    }
		
}