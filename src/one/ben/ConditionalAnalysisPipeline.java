package one.ben;

import gwas.FAST;
import gwas.FAST.DataDefinitions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import one.ScriptExecutor;

import common.Aliases;
import common.Array;
import common.Files;
import common.HashVec;
import common.Positions;
import common.ext;

public class ConditionalAnalysisPipeline {
    
    static class Region {
        String label;
        String indexSNP;
        int chr;
        int start;
        int stop;
        String analysisRootDir;
        String regionDirNameRoot;
        // set programmatically:
        String genoData;
        String infoData;
        
        HashSet<String> prevSNPs = new HashSet<String>();
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder("region[id:" + label + ", SNP:" + indexSNP + ", UCSC:chr" + chr + ":" + start + ":" + stop + "]");
            return sb.toString();    
        }
        
    }
    
    static class ConditionalAnalysisToolset_FAST implements Runnable {
        
        static final String DATA_FILE_TEMPLATE = "chr#.<1>.<2>";
        static final String DATA_FILE_DELIMITER = ".";
        static final String TEMP_GENO_FILE_EXTENSION = ".temp.gz";
        static final String TEMP_INFO_FILE_EXTENSION = ".temp_info";
        static final String MISSING_DATA = ".";
        static final int NUM_THREADS = 24;
        static final double PVAL_THRESHOLD = 0.0001;
        
        static final Object PRINT_LOCK = new Object();
        
        private static String[] findFiles(final Region region, final DataDefinitions dataDefs) {
            String[] chrDataFiles = (new File(dataDefs.dataDir)).list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    // TODO assuming datafile names start with "chr#."
                    return name.startsWith("chr" + region.chr + DATA_FILE_DELIMITER) && name.endsWith(dataDefs.dataSuffix);
                }
            });
            // TODO generic parsing for file-name template: is it in chr#.<>.<> format [chunked], or chr# format [whole_chr]?
            // assuming data file names include position
            ArrayList<String> inclDataFiles = new ArrayList<String>();
            for (int i = 0; i < chrDataFiles.length; i++) {
                String[] pts = chrDataFiles[i].split("\\.");
                int chunkStart = Integer.parseInt(pts[1]);
                int chunkStop = Integer.parseInt(pts[2]);
                boolean checkStartOverlap = region.start <= chunkStart && region.stop >= chunkStart;
                boolean checkOverlap1 = region.start >= chunkStart && region.stop <= chunkStop;
                boolean checkOverlap2 = region.start <= chunkStart && region.stop >= chunkStop;
                boolean checkStopOverlap = region.start < chunkStop && region.stop > chunkStop;
                if (checkStartOverlap || checkOverlap1 || checkOverlap2 || checkStopOverlap) {
                    inclDataFiles.add(chrDataFiles[i]); 
                }
            }
            synchronized(PRINT_LOCK) {
                System.out.println(ext.getTime() + "]\tFound " + inclDataFiles.size() + " out of " + chrDataFiles.length + " possibly-matching data files!");
            }
            return inclDataFiles.toArray(new String[inclDataFiles.size()]);
        }
        
        private static String extractGenoAndInfoDataForRegion(final Region region, final DataDefinitions dataDefs, String[] dataFiles, String tempDir) {
            // create subdir for study/pop
            // TODO currently DOES NOT reuse data files!  This majorly increases space used, but we have problems reading from shared files in a multi-threaded environment
            String tempDataDir = tempDir + region.label + "/" + region.indexSNP + "/" + dataDefs.study + "_" + dataDefs.popcode + "/";
            new File(tempDataDir).mkdirs();
            
            String infoHeader = "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";
            
            TreeSet<String> sortedFiles = new TreeSet<String>(new Comparator<String>() {
                @Override
                public int compare(String o1, String o2) {
                    String[] parts1 = o1.split("\\.");
                    String[] parts2 = o2.split("\\.");
                    int start1 = Integer.parseInt(parts1[1]);
//                    int stop1 = Integer.parseInt(parts1[2]);
                    int start2 = Integer.parseInt(parts2[1]);
//                    int stop2 = Integer.parseInt(parts2[2]);
                    return Integer.valueOf(start1).compareTo(Integer.valueOf(start2));
                }
            });
            for (String file : dataFiles) {
                sortedFiles.add(file);
            }
            
            String newGenoFileName = "chr" + region.chr + "." + region.start + "." + region.stop + ".temp.gz";
            String newInfoFileName = "chr" + region.chr + "." + region.start + "." + region.stop + ".temp_info";
            
            boolean found = false;
            if (Files.exists(tempDataDir + newGenoFileName) && Files.exists(tempDataDir + newInfoFileName)) {
                synchronized(PRINT_LOCK) {
                    System.out.println(ext.getTime() + "]\tWARNING - Skipping data file export, files already exist!");
                    System.out.println(ext.getTime() + "]\tLoading geno/info data for index SNP [" + region.indexSNP + "]...");
                }
                found = true;
            }
            
            PrintWriter genoWriter = null;
            PrintWriter infoWriter = null;
            if (!found) {
                genoWriter = Files.getAppropriateWriter(tempDataDir + newGenoFileName);
                infoWriter = Files.getAppropriateWriter(tempDataDir + newInfoFileName);
            
                infoWriter.println(infoHeader);
            }
            
            for (String file : sortedFiles) {
                try {
                    BufferedReader genoReader = Files.getAppropriateReader(dataDefs.dataDir + file);
                    BufferedReader infoReader = Files.getAppropriateReader(dataDefs.dataDir + file.substring(0, file.length() - 3) + "_info");
                    
                    String infoLine = infoReader.readLine();
                    String delim = ext.determineDelimiter(infoLine);
                    String genoLine = null;
                    int count = 0;
                    while((infoLine = infoReader.readLine()) != null && (genoLine = genoReader.readLine()) != null) {
                        count++;
                        // geno/info lines should be one to one
                        String[] infoParts = infoLine.split(delim);
                        if (infoParts[1].trim().equals(region.indexSNP)) {
                            region.genoData = genoLine;
                            region.infoData = infoLine;
                            if (found) {
                                break;
                            }
                        }
                        int mkrPos = Integer.parseInt(infoParts[2]);
                        if (mkrPos < region.start || mkrPos > region.stop) {
                            continue;
                        }
                        if (!found) {
                            genoWriter.println(genoLine);
                            infoWriter.println(infoLine);
                        }
                    }
                    System.out.println(ext.getTime() + "]\tRead " + count + " data/info lines...");
                    
                    genoReader.close();
                    infoReader.close();
                } catch (FileNotFoundException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
            if (!found) {
                infoWriter.flush();
                genoWriter.flush();
                infoWriter.close();
                genoWriter.close();
            }
            
            return tempDataDir;
        }
        
        private static void createNewTraitFiles(final Region region, String traitDir, DataDefinitions dd) {
            String[] iids = HashVec.loadFileToStringArray(dd.indivFile, false, new int[]{0}, false);
            HashMap<String, Integer> indexMap = new HashMap<String, Integer>();
            for (int i = 0; i < iids.length; i++) {
                indexMap.put(iids[i], i);
            }
            
            HashMap<String, HashMap<String, HashMap<String, String>>> studyToFactorToPopToFile = FAST.loadTraitFiles(traitDir);
            for (Entry<String, HashMap<String, HashMap<String, String>>> studyMap : studyToFactorToPopToFile.entrySet()) {
                String study = studyMap.getKey();
                if (!dd.study.equals(study)) {
                    continue;
                }
                for (Entry<String, HashMap<String, String>> factorMap : studyMap.getValue().entrySet()) {
                    for (Entry<String, String> popMap : factorMap.getValue().entrySet()) {
                        String pop = popMap.getKey();
                        if (!dd.popcode.equals(pop)) {
                            continue;
                        }
                        createNewTraitFile(region, traitDir, popMap.getValue(), indexMap);
                    }
                }
            }
            
        }
        
        private static void createNewTraitFile(Region region, String traitDir, String traitFile, HashMap<String, Integer> iids) {
            String[] pts = traitFile.substring(0, traitFile.lastIndexOf(".")).split("_");
            String study = pts[0];
            String pop = pts[1];
            String factor = pts[2];
            String newTraitFile = study + "_" + pop + "_" + factor/* + (ext.replaceWithLinuxSafeCharacters(region.indexSNP, false).replaceAll("_", ""))*/ + ".trait";
            String dir = region.analysisRootDir + region.regionDirNameRoot;
            
            int offset = 5; // column index offset to start of geno data
            String[] genoData = region.genoData.split("[\\s]+");
            
            ArrayList<String> missing = new ArrayList<String>();
            
            try {
                BufferedReader reader = Files.getAppropriateReader(traitDir + traitFile);
                PrintWriter writer = Files.getAppropriateWriter(dir + newTraitFile);
                
                String line = reader.readLine();
                line = line + "\t" + region.indexSNP;
                
                writer.println(line);
                
                while ((line = reader.readLine()) != null) {
                    String iid = line.split("\t")[1];
                    Integer iidIndex = iids.get(iid);
                    String geno = MISSING_DATA;
                    if (iidIndex == null) {
                        missing.add(iid);
                    } else {
                        int iidInd = iidIndex.intValue();
//                        double geno1 = Double.parseDouble(genoData[offset + (3 * iidInd)]);
                        double geno2 = Double.parseDouble(genoData[offset + (3 * iidInd) + 1]);
                        double geno3 = Double.parseDouble(genoData[offset + (3 * iidInd) + 2]);
                        geno = "" + (geno2 + (2 * geno3));
                    }
                    
                    writer.println(line + "\t" + geno);
                }
                writer.flush();
                writer.close();
                reader.close();
                
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tWARNING - " + missing.size() + " individuals from .trait file are missing genotype data."); }
                
            } catch (FileNotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        Region region;
        String dataFile;
        String traitDir;
        String tempDir;
        HashMap<String, HashMap<String, DataDefinitions>> dataDefs;
        
        public ConditionalAnalysisToolset_FAST(Region region, String dataFile, String traitDir, String tempDir, HashMap<String, HashMap<String, DataDefinitions>> dataDefs) {
            this.region = region;
            this.dataFile = dataFile;
            this.traitDir = traitDir;
            this.tempDir = tempDir;
            this.dataDefs = dataDefs;
        }
        
        @Override
        public void run() {
//            String[] dataFileAndTraitDir = new String[2];
            // Study -> PopCode -> Defs
            ArrayList<DataDefinitions> allDefs = new ArrayList<FAST.DataDefinitions>();
            for (HashMap<String, DataDefinitions> sub : dataDefs.values()) {
                allDefs.addAll(sub.values());
            }
            
            ArrayList<String> newDataDefs = new ArrayList<String>();
            for (DataDefinitions dd : allDefs) {
                
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tRetrieving required data files..."); }
                String[] files = findFiles(region, dd);
                
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tExtracting region-specific geno/info data..."); }
                String dir = extractGenoAndInfoDataForRegion(region, dd, files, tempDir);
                
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tCreating new .trait files..."); }
                createNewTraitFiles(region, traitDir, dd);
                
                StringBuilder newDef = new StringBuilder();
                newDef.append(dd.study).append("\t")
                        .append(dd.popcode).append("\t")
                        .append(dir).append("\t")
                        .append(".temp.gz").append("\t")
                        // TODO sex-specific
                        .append(dd.indivFile);
                newDataDefs.add(newDef.toString());
            }
            synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tWriting new data.txt file..."); }
            String newDataFile = region.analysisRootDir + region.regionDirNameRoot + "data_" + ext.replaceWithLinuxSafeCharacters(region.indexSNP, false) + ".txt";
            Files.writeList(newDataDefs.toArray(new String[newDataDefs.size()]), newDataFile);
            String regionPath = region.analysisRootDir + region.regionDirNameRoot;
            
            try {
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tPreparing FAST analysis in directory [" + regionPath + "]..."); }
                String[] analysisDirs = FAST.prepareFAST(regionPath, newDataFile, regionPath, false);
                
                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tRunning " + analysisDirs.length + " FAST analyses..."); }
                boolean[] runs = Array.booleanArray(analysisDirs.length, false);
                for (int i = 0; i < analysisDirs.length; i++) {
                    (new ScriptExecutor(NUM_THREADS)).run(analysisDirs[i] + "input.txt", "took");
                    runs[i] = ScriptExecutor.outLogExistsComplete(analysisDirs[i] + "output/input.log_0.out", "took");
                    if (!runs[i]) {
                        // TODO Error - FAST failed!
                    }
                }

                synchronized(PRINT_LOCK) { System.out.println(ext.getTime() + "]\tProcessing FAST results..."); }
                for (String study : dataDefs.keySet()) {
                    String studyDir = regionPath + study + "/";
                    
                    if (!Files.exists(studyDir)) {
                        System.err.println(ext.getTime() + "]\tERROR - analysis directory [" + studyDir + "] does not exist.");
                        continue;
                    }
                    
                    FAST.processAndPrepareMETAL(studyDir);
                    
                    HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);
                    
                    // should only ever be one directory...
                    String[] factorDirs = Files.listDirectories(studyDir, false);
                    for (String factorDir : factorDirs) {
                        String dir = studyDir + factorDir + "/";
                        String file = null;
                        if (popDefs.size() == 1) {
                            // one population, no meta analysis from which to pull results
                            file = popDefs.keySet().toArray(new String[1])[0] + "/output/concatenated.result"; 
                        } else {
                            file = factorDir + "_InvVar1.out";
                        }
                        
                        if (Files.exists(dir + file)) {
                            String newSNP = extractIndexSnp(dir + file, region, PVAL_THRESHOLD);
                            
                            if (newSNP == null) {
                                synchronized (PRINT_LOCK) {
                                    System.out.println(ext.getTime() + "]\tCouldn't find a candidate SNP for iterative analysis; recursive analysis for region [" + region.label + "] complete.");
                                }
                            } else {
                                synchronized (PRINT_LOCK) {
                                    System.out.println(ext.getTime() + "]\tIterating analysis with most-significant SNP [" + newSNP + "]");
                                }
                                Region r2 = new Region();
                                r2.chr = region.chr;
                                r2.start = region.start;
                                r2.stop = region.stop;
                                r2.label = region.label;
                                r2.indexSNP = newSNP;
    
                                String newDir = r2.label + "_iter" + (region.prevSNPs.size() + 1) + "_" + ext.replaceWithLinuxSafeCharacters(r2.indexSNP, false) + "/";
                                new File(region.analysisRootDir + newDir).mkdirs();
                                r2.analysisRootDir = region.analysisRootDir;
                                r2.regionDirNameRoot = newDir;
                                r2.prevSNPs.addAll(region.prevSNPs);
                                r2.prevSNPs.add(region.indexSNP);
                                (new ConditionalAnalysisToolset_FAST(r2, dataFile, regionPath, tempDir, dataDefs)).run();
                            }
                        } else {
                            synchronized(PRINT_LOCK) {
                                System.out.println(ext.getTime() + "]\tError - file [" + dir + file + "] not found!");
                            }
                            // TODO error message - result file not found!
                        }
                        
                    }
                    
                }
                
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        
        private static String extractIndexSnp(String resultFile, Region region, double thresh) throws IOException {
            String[][] aliases = {Aliases.MARKER_NAMES, Aliases.PVALUES};
            
            BufferedReader reader = Files.getAppropriateReader(resultFile);
            
            String line = reader.readLine();
            
            String delim = ext.determineDelimiter(line);
            String[] hdr = line.split(delim);
            int[] indices = ext.indexFactors(aliases, hdr, false, true, false, false);
            
            double minPVal = thresh;
//            String minPValStr = null;
            String minPValSNP = null;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split(delim);
                if (parts[indices[0]].equals(region.indexSNP) || region.prevSNPs.contains(parts[indices[0]])) {
                    continue;
                }
                double pval = thresh;
                try {
                    pval = Double.parseDouble(parts[indices[1]]);
                } catch (NumberFormatException e) {}
                if (pval < minPVal) {
                    minPVal = pval;
//                    minPValStr = parts[indices[1]];
                    minPValSNP = parts[indices[0]];
                }
            }
            reader.close();
            
            return minPValSNP;
        }
        
    }
    
    private Region[] parseSetupFile(String file) {
        String[][] rgnDefs = HashVec.loadFileToStringMatrix(file, false, new int[]{0, 1, 2}, false);
        Region[] rgns = new Region[rgnDefs.length];
        for (int i = 0; i < rgnDefs.length; i++) {
            Region rgn = new Region();
            rgn.label = rgnDefs[i][0];
            rgn.indexSNP = rgnDefs[i][1];
            int[] tmp = Positions.parseUCSClocation(rgnDefs[i][2]);
            rgn.chr = tmp[0];
            rgn.start = tmp[1];
            rgn.stop = tmp[2];
            rgns[i] = rgn;
        }
        return rgns;
    }
    
    private boolean[] buildAnalysisFolders(String dir, Region[] rgns) {
        boolean[] results = Array.booleanArray(rgns.length, false);
        for (int i = 0; i < rgns.length; i++) {
            String newDir = rgns[i].label + "_iter0_" + (ext.replaceWithLinuxSafeCharacters(rgns[i].indexSNP, false).replaceAll("_", ""))/* + "_" + rgns[i].chr + "_" + rgns[i].start + "_" + rgns[i].stop*/ + "/";
            results[i] = new File(dir + newDir).mkdirs();
            rgns[i].analysisRootDir = dir;
            rgns[i].regionDirNameRoot = newDir;
        }
        return results;
    }

    /*
     * Input file format:
     * Tab-delimited columns:
     * REGION_LABEL   INDEX_SNP   UCSC_REGION
     */
    private void setup(String analysisDir, String inputFile, String dataFile, String tempDataDir, String traitDir) {
        System.out.println(ext.getTime() + "]\tParsing regions from input file...");
        Region[] rgns = parseSetupFile(inputFile);
        System.out.println(ext.getTime() + "]\tFound " + rgns.length + " regions for analysis");
        boolean[] dirCreation = buildAnalysisFolders(analysisDir, rgns);
        System.out.println(ext.getTime() + "]\tParsing data file..."); // TODO divorce from FAST?
        HashMap<String, HashMap<String, DataDefinitions>> dataDefs = null;
        try {
            dataDefs = FAST.parseFile(dataFile);
        } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }

        int threads = rgns.length;
        int avail = Runtime.getRuntime().availableProcessors();
        threads = Math.min(threads, avail);
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        
        // TODO error checking on folder creation
        for (int i = 0; i < rgns.length; i++) {
            if (dirCreation[i]) {
                System.out.println(ext.getTime() + "]\tProcessing region " + rgns[i].toString());
                ConditionalAnalysisToolset_FAST run = new ConditionalAnalysisToolset_FAST(rgns[i], dataFile, traitDir, tempDataDir, dataDefs);
                executor.execute(run);
            }
        }

        executor.shutdown();
        try {
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        
        String analysisDir = "";
        String inputFile = "";
        String dataFile = "";
        String tempDataDir = "";
        String traitDir = "";
        
        String usage = "WRONG";
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("dir=")) {
                analysisDir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("input=")) {
                inputFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("data=")) {
                dataFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("temp=")) {
                tempDataDir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("trait=")) {
                traitDir = args[i].split("=")[1];
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.err.println(Array.toStr(args, "\n"));
            System.exit(1);
        }
        (new ConditionalAnalysisPipeline()).setup(analysisDir, inputFile, dataFile, tempDataDir, traitDir);
    }
    
    

//    String initialResultFile;
//    
//    String[][] RESULT_HEADER_COLS_LOOKUP = {Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS, Aliases.PVALUES};
//    String[][] HW_HEADER_COLS_LOOKUP = {Aliases.CHRS, {"RegionStart"}, {"RegionStop"}};
////    #Fam_ID Ind_ID  Dad_ID  Mom_ID  Sex     Phenotype       Age     Male    PC1     PC2     PC3     PC4     rs867186
//    
//    float DEFAULT_INDEX_THRESHOLD = 0.00000005f;
//    int DEFAULT_WIN_MIN_SIZE_PER_SIDE = 500000;
//    float DEFAULT_WIN_EXT_THRESHOLD = 0.000005f;
//    
//    PrintStream logStream = System.out;
////    logStream.println(ext.getTime() + "]\t");
//    
//    private void setup(String inputResultFile, String originalPhenoFile) {
//        this.initialResultFile = inputResultFile;
//        int[] cols = ext.indexFactors(RESULT_HEADER_COLS_LOOKUP, Files.getHeaderOfFile(initialResultFile, null), false, true, false, false);
//        
//        logStream.println(ext.getTime() + "]\tLoading marker data...");
//        String[][] resultFileData = HashVec.loadFileToStringMatrix(initialResultFile, true, cols, false);
//        
//        logStream.println(ext.getTime() + "]\tRunning HitWindows...");
//        String[][] windows = HitWindows.determine(initialResultFile, DEFAULT_INDEX_THRESHOLD, DEFAULT_WIN_MIN_SIZE_PER_SIDE, DEFAULT_WIN_EXT_THRESHOLD, new String[0]);
//        String[] hdr = windows[0];
//        int[] rgnIndices = ext.indexFactors(HW_HEADER_COLS_LOOKUP, hdr, false, true, false, false);
//        
//        logStream.println(ext.getTime() + "]\tFound " + (windows.length - 1) + " HitWindows; extracting min-p-value SNPs and constructing conditional sub-analyses...");
//        for (int i = 1; i < windows.length; i++) {
//            Window windowDtls = processHitWindow(windows[i], rgnIndices, resultFileData);
//        }
//        
//    }
//    
//    
//    private Window processHitWindow(String[] window, int[] indices, String[][] markerData) {
//        int chr = Integer.parseInt(window[indices[0]]);
//        int rgnStart = Integer.parseInt(window[indices[1]]);
//        int rgnStop = Integer.parseInt(window[indices[2]]);
//        
//        String[] minPValMarkerLine = null;
//        double minPval = Double.MAX_VALUE;
//        for (int i = 0; i < markerData.length; i++) {
////            String mkrName = markerData[i][0];
//            int mkrChr = Integer.parseInt(markerData[i][1]);
//            if (mkrChr < chr) continue;
//            if (mkrChr > chr) break;
//            int mkrPos = Integer.parseInt(markerData[i][2]);
//            if (mkrPos < rgnStart) continue;
//            if (mkrPos > rgnStop) break;
//            double mkrPval = Double.parseDouble(markerData[i][3]);
//            if (mkrPval < minPval) {
//                minPval = mkrPval;
//                minPValMarkerLine = markerData[i];
//            }
//        }
//        Window w = new Window();
//        w.chr = chr;
//        w.start = rgnStart;
//        w.stop = rgnStop;
//        w.minSNP = minPValMarkerLine;
//        return w;
//    }
    
}


