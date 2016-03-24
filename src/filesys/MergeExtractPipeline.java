package filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
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
import common.Positions;
import common.ext;

public class MergeExtractPipeline {
    
    private class DataSource {
        String label, dataFile, mapFile, idFile;
        boolean relD, relM, relI;
        public DataSource(String label, String dir, String d, String m, String i) {
            this.relD = Files.isRelativePath(d);
            this.relM = Files.isRelativePath(m);
            this.relI = Files.isRelativePath(i);
            this.dataFile = (relD ? dir : "") + d;
            this.mapFile = (relM ? dir : "") + m;
            this.idFile = (relI ? dir : "") + i;
            this.label = label;
        }
    }
    
    ArrayList<DataSource> dataSources;
    private String outFileD = null;
    private String outFileM = null;
    int outFormat = -1;
    String markersFile = null;
    String[] markers = null;
    int[][] markerLocations;
    String runDir = null;
    Logger log = null;
    String regionsFile = null;
    int[][] regions;
    String[] regionLabels;
    
    boolean splitOutput = false;
    boolean overwrite = false;
    boolean renameMarkers = true;
    
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
    
    /**
     * Either absolute path to file, or relative path to file in either previously-specified run-directory or current directory (depending on where the file exists)
     * @param regionFile
     * @return
     */
    public MergeExtractPipeline setRegions(String regionFile) {        
        BufferedReader reader;
        String line;
        String[] temp;
        boolean hasLabels = false;
        ArrayList<int[]> reg;
        ArrayList<String> lbl;
        
        initLog();

        if (this.markersFile != null) {
            String msg = "Error - cannot specify both markers-to-extract and regions-to-extract!";
            log.reportError(msg);
            throw new IllegalArgumentException(msg);
        }
        
        lbl = new ArrayList<String>();
        reg = new ArrayList<int[]>();
        this.regionsFile = regionFile;
        if (Files.isRelativePath(regionFile)) {
            if (runDir != null && !"".equals(runDir)) {
                if (Files.exists(runDir + regionFile)) {
                    this.regionsFile = runDir + regionFile;
                }
            } else {
                if (Files.exists("./" + regionFile)) {
                    this.regionsFile = "./" + regionFile;
                }
            }
        }
        try {
            reader = Files.getAppropriateReader(this.regionsFile);
            line = reader.readLine();
            if (line != null) {
                temp = line.split("\t");
                if (temp.length > 1) {
                    hasLabels = true;
                }
                reg.add(Positions.parseUCSClocation(temp[0]));
                if (hasLabels) {
                    lbl.add(temp[1]);
                } else {
                    lbl.add(ext.replaceWithLinuxSafeCharacters(temp[0], true));
                }
            }
            while ((line = reader.readLine()) != null) {
                temp = line.split("\t");
                reg.add(Positions.parseUCSClocation(temp[0]));
                if (hasLabels) {
                    lbl.add(temp[1]);
                } else {
                    lbl.add(ext.replaceWithLinuxSafeCharacters(temp[0], true));
                }
            }
            reader.close();
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
        
        this.regions = reg.toArray(new int[reg.size()][]);
        this.regionLabels = lbl.toArray(new String[lbl.size()]);
        
        log.report("Set regions-to-extract file: " + this.regionsFile);
        
        return this;
    }
    
    public MergeExtractPipeline setMarkers(String mkrsFile) {
        initLog();
        
        if (this.regionsFile != null) {
            String msg = "Error - cannot specify both markers-to-extract and regions-to-extract!";
            log.reportError(msg);
            throw new IllegalArgumentException(msg);
        }
        
        this.markersFile = mkrsFile;

        if (Files.isRelativePath(mkrsFile)) {
            if (runDir != null && !"".equals(runDir)) {
                if (Files.exists(runDir + mkrsFile)) {
                    this.markersFile = runDir + mkrsFile;
                }
            } else {
                if (Files.exists("./" + mkrsFile)) {
                    this.markersFile = "./" + mkrsFile;
                }
            }
        }
        this.markers = HashVec.loadFileToStringArray(this.markersFile, false, new int[]{0}, false);
        SnpMarkerSet markerSet = new SnpMarkerSet(markers);
        markerSet.parseSNPlocations(log);
        this.markerLocations = markerSet.getChrAndPositionsAsInts();
        for (int i = 0; i < markerLocations.length; i++) {
            if (markerLocations[i][0] == 0 || markerLocations[i][1] == -1) {
                log.reportError("Error - one or more markers could not be located.  Merging/Extracting may fail for this reason.");
                break;
            }
        }
        log.report("Set markers-to-extract file: " + this.markersFile);
        return this;
    }
    
    public MergeExtractPipeline setRenameMarkers(boolean rename) {
        this.renameMarkers = rename;
        return this;
    }

    public MergeExtractPipeline setRunDirectory(String runDir, boolean create) {
        initLog();
        this.runDir = runDir;
        if (!Files.exists(runDir)) {
            if (create) {
                (new File(runDir)).mkdirs();
            } else {
                throw new IllegalArgumentException("Error - specified run directory \"" + runDir + "\" doesn't exist, and create flag wasn't set.  Please fix or create directory and try again.");
            }
        }
        return this;
    }
    
    public MergeExtractPipeline setOutputFiles(String dataFile, String mapFile) {
        this.outFileD = dataFile;
        this.outFileM = mapFile;
        return this;
    }
    
    public MergeExtractPipeline setSplitOutputByRegions(boolean split) {
        this.splitOutput = split;
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

    public MergeExtractPipeline addDataSource(String lbl, String dir, String dataFile, String mapFile, String idFile) {
        initLog();
        dataSources.add(new DataSource(lbl, dir, dataFile, mapFile, idFile));
        log.report("Added data source: " + dataSources.get(dataSources.size() - 1).dataFile);
        return this;
    }
    
    public MergeExtractPipeline addDataSources(final String lbl, final String dir, final String dataFileExt, final String mapFileExt, final String idFile, final int bpWindow) {
        initLog();
        
        log.report("Looking for data files in directory \"" + dir + "\" that are named with the pattern \"chr##.##bpStart##.##bpEnd##.fileExt\".");
        log.reportTimeWarning("if data files are not named according to this scheme, all data files will be included.");

        String[] filesToAdd = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                boolean keep = false;
                if (name.endsWith(dataFileExt)) {
                    if (/*markerLocations == null && */regions == null) {
                        keep = true;
                    } else {
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
                            if (markerLocations != null) {
                                for (int i = 0; i < markerLocations.length; i++) {
                                    if (markerLocations[i][0] == chr) {
                                        if (bpStart != -1 && bpEnd != -1) {
                                            if (bpStart <= markerLocations[i][1] - bpWindow && bpEnd >= markerLocations[i][1] + bpWindow) { 
                                                keep = true;
                                                break;
                                            }
                                        } else {
                                            keep = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            if (regions != null) {
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
            addDataSource(lbl, dir, file, file.substring(0, file.length() - dataFileExt.length()) + mapFileExt, idFile);
            log.report("Added data source: " + dataSources.get(dataSources.size() - 1).dataFile);
        }
        return this;
    }
    
    private String getOutputDataFile() {
        return (Files.isRelativePath(outFileD) ? (runDir == null ? "./" : runDir) : "") + outFileD;
    }

    private String getOutputDataFile(String subDir) {
        return ext.verifyDirFormat((Files.isRelativePath(outFileD) ? (runDir == null ? "./" : runDir) : "") + subDir) + outFileD;
    }

    private String getOutputMapFile() {
        return (Files.isRelativePath(outFileM) ? (runDir == null ? "./" : runDir) : "") + outFileM;
    }

    private String getOutputMapFile(String subDir) {
        return ext.verifyDirFormat((Files.isRelativePath(outFileM) ? (runDir == null ? "./" : runDir) : "") + subDir) + outFileM;
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
        if (markersFile == null) {
            log.reportTimeWarning("no markers file specified; all markers will be exported.");
        } else if (!Files.exists(markersFile)) {
            log.reportTimeError("specified markers file not found; please check file path and try again.");
            return;
        }

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
        ArrayList<String> plinkLabels = new ArrayList<String>();
        // discover and merge all plink files  
        for (int i = dataSources.size() - 1; i >= 0; i--) {
            if (dataSources.get(i).dataFile.endsWith(".bed") || dataSources.get(i).dataFile.endsWith(".ped")) {
                plinkRoots.add(ext.rootOf(dataSources.get(i).dataFile, false));
                plinkLabels.add(dataSources.get(i).label);
                dataSources.remove(i);
            }
        }
        if (plinkRoots.size() > 0) {
            String[] roots = Array.toStringArray(plinkRoots);
            String[] lbls = Array.toStringArray(plinkLabels);
            plinkRoots = null;
            String outRoot = dir + "plink_merged";
            String mergeCommand = PlinkMergePrep.merge(PlinkMergePrep.BEDBIMFAM, outRoot, overwrite, renameMarkers, regionsFile, markersFile, roots, lbls);

            log.report("Running PLINK merge command: " + mergeCommand);
            boolean result = CmdLine.runDefaults(mergeCommand, dir);
            if (!result) {
                log.reportTimeError("PLINK merge command failed.  Please check output for errors and try again.");
                return;
            }
            dataSources.add(new DataSource(null, dir, PSF.Plink.getBED(outRoot), PSF.Plink.getBIM(outRoot), PSF.Plink.getFAM(outRoot))); // no prepend here, as we've already renamed the markers using PlinkMergePrep
            if (markersFile != null && !"".equals(markersFile)) {
                String newMkrFile = (new File(PlinkMergePrep.TEMP_MKR_FILE)).getAbsolutePath();
                log.report("Setting markers file to temporarily generated plink-renamed markers file: " + newMkrFile);
                this.markersFile = newMkrFile;
            }
        }
        
        System.gc();
        
        log.reportTime("Merging data...");
        log.reportTime("Starting from data file: " + dataSources.get(0).dataFile);
        DosageData dd1 = new DosageData(dataSources.get(0).dataFile, dataSources.get(0).idFile, dataSources.get(0).mapFile, regionsFile, markersFile, renameMarkers ? dataSources.get(0).label : "", true, null);
        for (int i = 1; i < dataSources.size(); i++) {
            System.gc();
            log.reportTime("... merging with data file: " + dataSources.get(i).dataFile);
            DosageData dd2 = new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile, dataSources.get(i).mapFile, regionsFile, markersFile, renameMarkers ? dataSources.get(i).label : "", true, null);
            dd1 = DosageData.combine(dd1, dd2, log);
            dd2 = null;
            System.gc();
        }
        
        log.reportTime("Writing final merged data to output file: " + getOutputDataFile());
        
        if (outFormat == -1) {
            outFormat = DosageData.determineType(getOutputDataFile());
        }
        if (splitOutput) {
            String outD, outM;
            for (int i = 0; i < regions.length; i++) {
                outD = getOutputDataFile(regionLabels[i]);
                outM = getOutputMapFile(regionLabels[i]);
                (new File(outD)).mkdirs();
                (new File(outM)).mkdirs();
                dd1.writeToFile(outD, outM, markers, new int[][]{regions[i]}, true, true, DosageData.PARAMETERS[outFormat], log);
            }
        } else {
            dd1.writeToFile(getOutputDataFile(), getOutputMapFile(), markers, regions, true, true, DosageData.PARAMETERS[outFormat], log);
        }
        
        System.gc();
    }
    
    private boolean checkMemory() {
        initLog();
        log.reportError("Warning - memory check not implemented yet - if not enough memory is provided, an Out Of Memory exception may occur.");
        return true;
    }
    
    private static void parseDataFile(MergeExtractPipeline mep, String data) {
        BufferedReader reader; 
        String line;
        String[] temp;
        
        if (null == data || "".equals(data) || !Files.exists(data)) {
            throw new IllegalArgumentException("Error - provided data file \"" + data + "\" doesn't exist.");
        }
        
        try {
            reader = Files.getAppropriateReader(Files.exists(mep.runDir + data) ? mep.runDir + data : "./" + data);
            while ((line = reader.readLine()) != null) {
                temp = line.split("\t");
                if (temp.length == 4) {
                    mep.addDataSource(temp[0], null, temp[1], temp[2], temp[3]);
                } else if (temp.length == 5) {
                    mep.addDataSources(temp[0], temp[1], temp[2], temp[3], temp[4], 0);
                } else {
                    System.err.println("Error - skipping invalid entry in data file: " + line);
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(e);
        }
    }

    //    MergeExtractPipeline pipeline = new MergeExtractPipeline();
//    pipeline.setMarkers(markersFile);           
//    pipeline.setRunDirectory("/scratch.global/cole0482/merge/", true);
//    pipeline.setOutputFormat(DosageData.DATABASE_DOSE_FORMAT);
//    pipeline.setOutputFiles(outFile, mapOutFile);
//    pipeline.setRenamePlinkMarkers(true);
//    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "gwas.bed", "gwas.bim", "gwas.fam");
//    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "exome.bed", "exome.bim", "exome.fam");
//    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "metab.bed", "metab.bim", "metab.fam");
    // add more data;
//    pipeline.run();
    public static void main(String[] args) {
        int numArgs = args.length;
        String outfileD = null;
        String outfileM = null;
        String data = "./dataSources.txt";
        String markers = null;
        String regions = null;
        String rundir = "./";
        String logFile = null;
        boolean split = false;
        boolean overwrite = false;
        boolean rename = true;
        
        String usage = "\n" + 
                       "filesys.MergeExtractPipeline requires 4+ arguments\n" + 
                       "   (1) Output Data filename (i.e. outData=" + outfileD + " (default))\n" + 
                       "   (2) Output Map filename (i.e. outMap=" + outfileM + " (default))\n" +
                       "   (3) Run directory (output files and temporary files will be created here) (i.e. runDir=" + rundir + " (default))\n" +
                       "   (4) File listing data sources (i.e. data=" + data + " (default))\n" + 
                       "          Example:\n" +
                       "          dataLabel1\tfullPathDataFile1\tFullPathMapFile1\tFullPathIdFile1\n" + 
                       "          dataLabel2\tfullPathDataFile2\tFullPathMapFile2\tFullPathIdFile2\n" + 
                       "          dataLabel3\tdir1\tdataFileExt1\tmapFileExt1\tidFile3\n" + 
                       "          dataLabel4\tdir2\tdataFileExt2\tmapFileExt2\tidFile4\n" + 
                       "   (5a) Regions-to-extract filename (i.e. regions=" + regions + " (default))\n" + 
                       "   (5b) Markers-to-extract filename (i.e. markers=" + markers + " (default))\n" +
                       "          (Note: only one is allowed, either regions or markers, not both)\n" + 
                       "   (6) Optional: Log file name (i.e. log=" + logFile + " (default))\n" +
                       "   (7) Optional: Split output files (if region file is specified) (i.e. split=" + split + " (default))\n" +
                       "   (8) Optional: Overwrite files if they already exist (i.e. overwrite=" + overwrite + " (default))\n" +
                       "   (9) Optional: Rename markers in any data files to dataLabel+MarkerName (i.e. rename=" + rename + " (default))\n" +
                       "\n";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("outData=")) {
                outfileD = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("outMap=")) {
                outfileM = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("runDir=")) {
                rundir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("markers=")) {
                markers = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("data=")) {
                data = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("regions=")) {
                regions = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("split=")) {
                split = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith("overwrite=")) {
                overwrite = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith("rename=")) {
                rename = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                logFile = args[i].split("=")[1];
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0 || args.length == 0) {
            System.err.println(usage);
            System.exit(1);
        }
        
        MergeExtractPipeline mep = new MergeExtractPipeline();
        if (logFile != null) { mep.setLogger(new Logger(logFile)); }
        if (rundir != null) { mep.setRunDirectory(rundir, true); }
        if (regions != null) { mep.setRegions(regions); }
        if (markers != null) { mep.setMarkers(markers); }
        if (overwrite) { mep.setOverwrite(); }
        mep.setSplitOutputByRegions(split);
        mep.setRenameMarkers(rename);
        mep.setOutputFiles(outfileD, outfileM);
        
        parseDataFile(mep, data);
        
        mep.run();
    }
    
}
