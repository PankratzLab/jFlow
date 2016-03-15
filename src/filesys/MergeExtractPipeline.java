package filesys;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import cnv.manage.PlinkMergePrep;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.ext;

public class MergeExtractPipeline {
    
    private class DataSource {
        String dataFile, mapFile, idFile;
        boolean relD, relM, relI;
        public DataSource(String dir, String d, String m, String i) {
            this.relD = Files.isRelativePath(d);
            this.relM = Files.isRelativePath(m);
            this.relI = Files.isRelativePath(i);
            this.dataFile = (relD ? dir : "") + d;
            this.mapFile = (relM ? dir : "") + m;
            this.idFile = (relI ? dir : "") + i;
        }
    }
    
    ArrayList<DataSource> dataSources;
    private String outFileD = null;
    private String outFileM = null;
    int outFormat = -1;
//    String markersFile = null;
//    String[] markers = null;
//    int[][] markerLocations;
    String runDir = null;
    Logger log = null;
    String regionsFile = null;
    int[][] regions;
    
    boolean overwrite = false;
    boolean renamePlink = true;
    
    public MergeExtractPipeline() {
        dataSources = new ArrayList<MergeExtractPipeline.DataSource>();
    }
    
    private void initLog() {
        if (log == null) {
            log = new Logger();
        }
    }
    
    public MergeExtractPipeline setLogger(Logger log) {
        this.log = log;
        return this;
    }
    
//    public MergeExtractPipeline setMarkers(String mkrsFile) {
//        initLog();
//        this.markersFile = mkrsFile;
//        this.markers = HashVec.loadFileToStringArray(mkrsFile, false, new int[]{0}, false);
//        SnpMarkerSet markerSet = new SnpMarkerSet(markers);
//        markerSet.parseSNPlocations(log);
//        this.markerLocations = markerSet.getChrAndPositionsAsInts();
//        for (int i = 0; i < markerLocations.length; i++) {
//            if (markerLocations[i][0] == 0 || markerLocations[i][1] == -1) {
//                log.reportError("Error - one or more markers could not be located.  Merging/Extracting may fail for this reason.");
//                break;
//            }
//        }
//        log.report("Set markers-to-extract file: " + this.markersFile);
//        return this;
//    }
    
    public MergeExtractPipeline setRenamePlinkMarkers(boolean rename) {
        this.renamePlink = rename;
        return this;
    }

    public MergeExtractPipeline setRunDirectory(String runDir) {
        this.runDir = runDir;
        return this;
    }
    
    /**
     * 
     * @param dataFile
     * @param mapFile
     * @return
     */
    public MergeExtractPipeline setOutputFiles(String dataFile, String mapFile) {
        this.outFileD = dataFile;
        this.outFileM = mapFile;
        return this;
    }

    /**
     * Refers to DosageData formats
     * 
     * @param format
     * @return
     */
    public MergeExtractPipeline setOutputFormat(int format) {
        this.outFormat = format;
        return this;
    }

    public MergeExtractPipeline setOverwrite() {
        this.overwrite = true;
        return this;
    }

    public MergeExtractPipeline addDataSource(String dir, String dataFile, String mapFile, String idFile) {
        initLog();
        dataSources.add(new DataSource(dir, dataFile, mapFile, idFile));
        log.report("Added data source: " + dataSources.get(dataSources.size() - 1).dataFile);
        return this;
    }
    
    public MergeExtractPipeline addDataSources(final String dir, final String dataFileExt, final String mapFileExt, final String idFile, final int bpWindow) {
        initLog();
        String[] filesToAdd = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                boolean keep = false;
                if (name.endsWith(dataFileExt)) {
//                    if (markerLocations == null) {
                    if (regions == null) {
                        keep = true;
                    } else {
                        log.report("Looking for data files in directory that are named with the pattern \"chr#.##bpStart##.##bpEnd##.fileExtension\".");
                        log.reportTimeWarning("if data files are not named according to this scheme, all data files will be included.");
                        Matcher m = Pattern.compile(DosageData.CHR_REGEX).matcher(name);
                        if (m.matches()) {
                            byte chr = -1;
                            chr = (byte) Integer.parseInt(m.group(1));
                            int bpStart = -1, bpEnd = -1;
                            String[] testPos = name.split("\\.");
                            for (int i = 0; i < testPos.length; i++) {
                                if (testPos[i].startsWith("chr")) {
                                    if (testPos.length > (i+2)) {
                                        bpStart = Integer.parseInt(testPos[i + 1]);
                                        bpEnd = Integer.parseInt(testPos[i + 2]);
                                    } else {
                                        log.reportError("not enough information in file name to determine bpStart/bpEnd from filename \"" + name + "\" - whole chromosomes will be included in data merge.");
                                    }
                                    break;
                                }
                            }
//                            for (int i = 0; i < markerLocations.length; i++) {
//                                if (markerLocations[i][0] == chr) {
//                                    if (bpStart != -1 && bpEnd != -1) {
//                                        if (bpStart <= markerLocations[i][1] - bpWindow && bpEnd >= markerLocations[i][1] + bpWindow) { 
//                                            keep = true;
//                                            break;
//                                        }
//                                    } else {
//                                        keep = true;
//                                        break;
//                                    }
//                                }
//                            }
                            for (int i = 0; i < regions.length; i++) {
                                if (regions[i][0] == chr) {
                                    if (bpStart != -1 && bpEnd != -1) {
                                        if (bpStart <= regions[i][1] - bpWindow && bpEnd >= regions[i][1] + bpWindow) { 
                                            keep = true;
                                            break;
                                        }
                                    } else {
                                        keep = true;
                                        break;
                                    }
                                }
                            }
                        } else {
                            keep = true;
                        }
                    }
                }
                return keep;
            }
        });
        for (String file : filesToAdd) {
            addDataSource(dir, file, file.substring(0, file.length() - dataFileExt.length()) + mapFileExt, idFile);
            log.report("Added data source: " + dataSources.get(dataSources.size() - 1).dataFile);
        }
        return this;
    }
    
    private String getOutputDataFile() {
        return (Files.isRelativePath(outFileD) ? (runDir == null ? "./" : runDir) : "") + outFileD;
    }

    private String getOutputMapFile() {
        return (Files.isRelativePath(outFileM) ? (runDir == null ? "./" : runDir) : "") + outFileM;
    }
    
    public void run() {
        initLog();
        if (dataSources.size() == 0) {
            log.reportTimeError("no data sources specified!");
            return;
        }
        if (getOutputDataFile() == null || "".equals(getOutputDataFile())) {
            log.reportTimeError("output data file name not specified!");
            return;
        }
        if (Files.exists(getOutputDataFile()) && !overwrite) {
            log.reportTimeError("output data file already exists!  Please remove file or set -overwrite flag, and re-run.");
            return;
        }
        if (getOutputMapFile() == null || "".equals(getOutputMapFile())) {
            log.reportTimeError("output map file name not specified!");
            return;
        }
        if (Files.exists(getOutputMapFile()) && !overwrite) {
            log.reportTimeError("output map file already exists!  Please remove file or set -overwrite flag, and re-run.");
            return;
        }
        if (regionsFile == null) {
            log.reportTimeWarning("no regions file specified; all data will be exported.");
        } else if (!Files.exists(regionsFile)) {
            log.reportTimeError("specified regions file not found; please check file path and try again.");
            return;
        }
//        if (markersFile == null) {
//            log.reportTimeWarning("no markers file specified; all markers will be exported.");
//        } else if (!Files.exists(markersFile)) {
//            log.reportTimeError("specified markers file not found; please check file path and try again.");
//            return;
//        }

        if (!checkMemory()) {
            return;
        }
        
        String dir = "./";
        if (runDir == null) {
            log.reportTimeWarning("no run directory specified, merging will take place in the current directory.");
        } else {
            dir = runDir;
        }
        dir = ext.verifyDirFormat((new File(dir)).getAbsolutePath());
        
//        // all relative paths are already absolute
//        for (int i = 0; i < dataSources.size(); i++) {
//            if (Files.isRelativePath(dataSources.get(i).dataFile)) {
//                dataSources.get(i).dataFile = dir + dataSources.get(i).dataFile;
//            }
//            if (Files.isRelativePath(dataSources.get(i).mapFile)) {
//                dataSources.get(i).mapFile = dir + dataSources.get(i).mapFile;
//            }
//            if (Files.isRelativePath(dataSources.get(i).idFile)) {
//                dataSources.get(i).idFile = dir + dataSources.get(i).idFile;
//            }
//        }
        
        ArrayList<String> plinkRoots = new ArrayList<String>();
        // discover and merge all plink files  
        for (int i = dataSources.size() - 1; i >= 0; i--) {
            if (dataSources.get(i).dataFile.endsWith(".bed") || dataSources.get(i).dataFile.endsWith(".ped")) {
                plinkRoots.add(ext.rootOf(dataSources.get(i).dataFile, false));
                dataSources.remove(i);
            }
        }
        if (plinkRoots.size() > 0) {
            String[] roots = Array.toStringArray(plinkRoots);
            plinkRoots = null;
            String outRoot = dir + "plink_merged";
            String mergeCommand = PlinkMergePrep.merge(PlinkMergePrep.BEDBIMFAM, outRoot, overwrite, renamePlink, regionsFile, null/*markersFile*/, roots);

            log.report("Running PLINK merge command: " + mergeCommand);
            boolean result = CmdLine.runDefaults(mergeCommand, dir);
            if (!result) {
                log.reportTimeError("PLINK merge command failed.  Please check output for errors and try again.");
                return;
            }
            dataSources.add(new DataSource(dir, PSF.Plink.getBED(outRoot), PSF.Plink.getBIM(outRoot), PSF.Plink.getFAM(outRoot)));
//            String newMkrFile = (new File(PlinkMergePrep.TEMP_MKR_FILE)).getAbsolutePath();
//            log.report("Setting markers file to temporarily generated plink-renamed markers file: " + newMkrFile);
//            this.markersFile = newMkrFile;
        }
        
        System.gc();
        
        log.reportTime("Merging data...");
        log.reportTime("Starting from data file: " + dataSources.get(0).dataFile);
        DosageData dd1 = new DosageData(dataSources.get(0).dataFile, dataSources.get(0).idFile, dataSources.get(0).mapFile, regionsFile, null/*markersFile*/, true, null);
        for (int i = 1; i < dataSources.size(); i++) {
            log.reportTime("... merging with data file: " + dataSources.get(i).dataFile);
            DosageData dd2 = new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile, dataSources.get(i).mapFile, regionsFile, null/*markersFile*/, true, null);
            dd1 = DosageData.combine(dd1, dd2);
            dd2 = null;
            System.gc();
        }
        
        log.reportTime("Writing final merged data to output file: " + getOutputDataFile());
        if (outFormat != -1) {
            dd1.writeToFile(getOutputDataFile(), getOutputMapFile(), null, true, DosageData.PARAMETERS[outFormat], log);
        } else {
            dd1.writeToFile(getOutputDataFile(), getOutputMapFile(), null, true, DosageData.PARAMETERS[DosageData.determineType(getOutputDataFile())], log);
        }
    }
    
    private boolean checkMemory() {
        initLog();
        log.reportError("Warning - memory check not implemented yet - if not enough memory is provided, an Out Of Memory exception may occur.");
        return true;
    }
    
    /*  
          String tgtMkrs = "";
          String runDir = "";
          String outData = "";
          String outMap = "";
          
          MergeExtractPipeline mep = new MergeExtractPipeline()
              .setMarkers(tgtMkrs)
              .setRunDirectory(runDir) // OPTIONAL
              .setOutputFiles(outData, outMap)
              .setOutputFormat(outFormat) // OPTIONAL
              .setOverwrite(true) // OPTIONAL
              .addDataSource()
              .addDataSource()
              .addDataSource()
              .addDataSources()
              .run();
              
     */
    
    
    
    
    
    
    
}
