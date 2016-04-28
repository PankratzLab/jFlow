package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;

import cnv.manage.ExportCNVsToPedFormat;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

public class EmimPipeline {
    
    static class PopFileData {
        String[] pops;
        LinkedHashMap<String, Integer> idsIndexMap;
        boolean[][] inclArr;
    }
    
    private static PopFileData loadPopFile(String popFile) throws IOException {
        BufferedReader reader = Files.getAppropriateReader(popFile);
        String line = reader.readLine();
        String[] hdrs = line.split("\t", -1);
        ArrayList<String> pops = new ArrayList<String>();
        ArrayList<Integer> popIndices = new ArrayList<Integer>();
        for (int i = 0; i < hdrs.length; i++) {
            if (!"".equals(hdrs[i])) {
                popIndices.add(i);
                pops.add(hdrs[i]);
            }
        }
        ArrayList<String> ids = new ArrayList<String>();
        @SuppressWarnings("unchecked")
        ArrayList<Boolean>[] incls = new ArrayList[pops.size()];
        for (int i = 0; i < incls.length; i++) {
            incls[i] = new ArrayList<Boolean>();
        }
        while ((line = reader.readLine()) != null) {
            String[] incl = line.split("\t", -1);
            for (int i = 0; i < incls.length; i++) {
                incls[i].add(Integer.parseInt(incl[popIndices.get(i)]) == 1);
            }
            ids.add(incl[0] + "\t" + incl[1]);
        }
        reader.close();
        reader = null;
        
        PopFileData popFileObj = new PopFileData();
        popFileObj.pops = pops.toArray(new String[pops.size()]);
        popFileObj.idsIndexMap = new LinkedHashMap<String, Integer>();
        for (int i = 0; i < ids.size(); i++) {
            popFileObj.idsIndexMap.put(ids.get(i), i);
        }
        popFileObj.inclArr = new boolean[incls.length][pops.size()];
        for (int i = 0; i < incls.length; i++) {
            popFileObj.inclArr[i] = Array.toBooleanArray(incls[i]);
        }
        
        hdrs = null;
        ids = null;
        incls = null;
        pops = null;
        popIndices = null;
        
        return popFileObj;
    }

    private static void generateKeepsFile(String keepsFile, PopFileData popData, int pop, PopFileData subPopData, int subPop) {
        PrintWriter writer = Files.getAppropriateWriter(keepsFile);
        for (String s : popData.idsIndexMap.keySet()) {
            if (popData.inclArr[pop][popData.idsIndexMap.get(s)]) {
                if (subPopData == null || (subPopData.idsIndexMap.containsKey(s) && subPopData.inclArr[subPop][subPopData.idsIndexMap.get(s)])) {
                    writer.println(s);
                }
            }
        }
        writer.flush();
        writer.close();
    }
    
    
    private static void generateFolderStructureAndKeepsFiles(String runDir, String[] cnvFiles, String[] plinkRoots, PopFileData popFile, PopFileData subPopFile) {
        if (cnvFiles != null) {
            for (String cnvFile : cnvFiles) {
                String cnvRoot = ext.rootOf(cnvFile, true);
                String cnvDir = runDir + cnvRoot + "/";
                boolean created = (!Files.exists(cnvDir) && (new File(cnvDir)).mkdir());
                if (!created) { 
                    /* TODO ERROR, or already exists */ 
                    System.err.println("Error - Could not create folder " + cnvDir);
                    continue; 
                }
                
                for (int p = 0; p < popFile.pops.length; p++) {
                    String popDir = cnvDir + ext.replaceWithLinuxSafeCharacters(popFile.pops[p], true) + "/";
                    created = (!Files.exists(popDir) && (new File(popDir)).mkdir());
                    if (!created) { 
                        /* TODO ERROR, or already exists */ 
                        System.err.println("Error - Could not create folder " + popDir);
                        continue; 
                    }
                    String keepsFile = popDir + "keeps.txt";
                    generateKeepsFile(keepsFile, popFile, p, null, 0);
                    
                    for (int sP = 0; sP < subPopFile.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopFile.pops[sP], true) + "/";
                        created = (!Files.exists(subPopDir) && (new File(subPopDir)).mkdir());
                        if (!created) { 
                            /* TODO ERROR, or already exists */ 
                            System.err.println("Error - Could not create folder " + subPopDir);
                            continue; 
                        }
                        keepsFile = subPopDir + "keeps.txt";
                        generateKeepsFile(keepsFile, popFile, p, subPopFile, sP);
                    }
                }
            }
        }
        if (plinkRoots != null) {
            for (String plinkRoot : plinkRoots) {
                String plinkDir = runDir + plinkRoot + "/";
                boolean created = (!Files.exists(plinkDir) && (new File(plinkDir)).mkdir());
                if (!created) { 
                    /* TODO ERROR, or already exists */ 
                    System.err.println("Error - Could not create folder " + plinkDir);
                    continue; 
                }
                
                for (int p = 0; p < popFile.pops.length; p++) {
                    String popDir = plinkDir + ext.replaceWithLinuxSafeCharacters(popFile.pops[p], true) + "/";
                    created = (!Files.exists(popDir) && (new File(popDir)).mkdir());
                    if (!created) { 
                        /* TODO ERROR, or already exists */ 
                        System.err.println("Error - Could not create folder " + popDir);
                        continue; 
                    }
                    String keepsFile = popDir + "keeps.txt";
                    generateKeepsFile(keepsFile, popFile, p, null, 0);
                    
                    for (int sP = 0; sP < subPopFile.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopFile.pops[sP], true) + "/";
                        created = (!Files.exists(subPopDir) && (new File(subPopDir)).mkdir());
                        if (!created) { 
                            /* TODO ERROR, or already exists */ 
                            System.err.println("Error - Could not create folder " + subPopDir);
                            continue; 
                        }
                        keepsFile = subPopDir + "keeps.txt";
                        generateKeepsFile(keepsFile, popFile, p, subPopFile, sP);
                    }
                }
            }
        }
    }
    
    static void setup(String runDir, String[] cnvFiles, String[] plinkRoots, String pedFile, String popFile, String subPopFile, double pThreshold, String qsubQueue, Logger log1) {
        ArrayList<String> pbsFiles = new ArrayList<String>();
        Logger log = log1 == null ? new Logger() : log1;
        PopFileData popData, subPopData;
        try {
            popData = loadPopFile(popFile);
            subPopData = loadPopFile(subPopFile);
        } catch (IOException e) {
            e.printStackTrace();
            log.reportException(e);
            return;
        }
        log.reportTimeInfo("Generating folder structure in run directory...");
        generateFolderStructureAndKeepsFiles(runDir, cnvFiles, plinkRoots, popData, subPopData);
        if (!checkRoots(cnvFiles, plinkRoots)) {
            log.reportError("Error - cannot have cnv roots and plink roots that have the same name.");
            return;
        }
        if (cnvFiles != null) {
            for (String cnvFile : cnvFiles) {
                String cnvRoot = ext.rootOf(cnvFile, true);
                String cnvDir = cnvRoot + "/";
                String plinkRoot = runDir + cnvRoot + "_0";
                if (!Files.exists(plinkRoot + ".bim") || !Files.exists(plinkRoot + ".bed") || !Files.exists(plinkRoot + ".fam")) {
                    log.reportTimeInfo("Exporting cnv files for root " + plinkRoot + " to PLINK files...");
                    ExportCNVsToPedFormat.export(cnvFile, pedFile, runDir + cnvRoot, "\r\n", ExportCNVsToPedFormat.PLINK_BINARY_FORMAT, true, true, false, false, false, false, Integer.MAX_VALUE, 0, log);
                }
                if (!Files.exists(plinkRoot + ".bim") || !Files.exists(plinkRoot + ".bed") || !Files.exists(plinkRoot + ".fam")) {
                    log.reportError("ERROR - couldn't find exported PLINK files for CNV root " + cnvRoot + " in directory " + runDir);
                    continue;
                }
                String relativePlinkRoot = "../" + cnvRoot + "_0";
                String resultFile = cnvRoot;
                
                log.reportTimeInfo("Generating EMIM files for cnv root " + plinkRoot);
                Emim.scriptAllInDir(runDir + cnvDir, plinkRoot, relativePlinkRoot, "GEN", null, pThreshold, resultFile);
                String pbsFile = cnvDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                pbsFiles.add(pbsFile);
    
                for (int p = 0; p < popData.pops.length; p++) {
                    String popDir = cnvDir + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "/";
                    resultFile = cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true);
                    Emim.scriptAllInDir(runDir + popDir, plinkRoot, "../" + relativePlinkRoot, "GEN", "keeps.txt", pThreshold, resultFile);
                    pbsFile = popDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                    if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                    pbsFiles.add(pbsFile);
                    
                    for (int sP = 0; sP < subPopData.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "/"; 
                        resultFile = cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_" + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true);
                        Emim.scriptAllInDir(runDir + subPopDir, plinkRoot, "../../" + relativePlinkRoot, "GEN", "keeps.txt", pThreshold, resultFile);
                        pbsFile = subPopDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                        if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                        pbsFiles.add(pbsFile);
                    }
                }
            }
        }
        if (plinkRoots != null) {
            for (String plinkRoot : plinkRoots) {
                String plinkDir = plinkRoot + "/";
                String resultFile = plinkRoot;

                log.reportTimeInfo("Generating EMIM files for PLINK root " + plinkRoot);
                Emim.scriptAllInDir(runDir + plinkDir, plinkRoot, "../" + plinkRoot, "GEN", null, pThreshold, resultFile);
                String pbsFile = plinkDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                pbsFiles.add(pbsFile);
                
                for (int p = 0; p < popData.pops.length; p++) {
                    String popDir = plinkDir + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "/";
                    resultFile = plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true);
                    Emim.scriptAllInDir(runDir + popDir, plinkRoot, "../../" + plinkRoot, "GEN", "keeps.txt", pThreshold, resultFile);
                    pbsFile = popDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                    if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                    pbsFiles.add(pbsFile);
                    
                    for (int sP = 0; sP < subPopData.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "/"; 
                        resultFile = plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_" + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true);
                        Emim.scriptAllInDir(runDir + subPopDir, plinkRoot, "../../../" + plinkRoot, "GEN", "keeps.txt", pThreshold, resultFile);
                        pbsFile = subPopDir + ext.rootOf(plinkRoot, true) + "_runEmim.pbs";
                        if (!Files.exists(runDir + pbsFile)) { /* TODO ERROR */ continue; }
                        pbsFiles.add(pbsFile);
                    }
                }
            }
        }
        

        log.reportTimeInfo("Writing QSUB files...");
        if (pbsFiles.size() > 0) {
            writeQsubs(pbsFiles, runDir, qsubQueue);
//            for (String pbsFile : pbsFiles) {
//                CmdLine.run("qsub " + (qsubQueue != null ? "-q " + qsubQueue : "") + pbsFile, ext.parseDirectoryOfFile(pbsFile));
//            }
            
        }
        
        String processCommand = "jcp gwas.EmimPipeline -process pThreshold=0.00001 dir=" + runDir;
        if (cnvFiles != null) {
            processCommand += " cnvs=" + Array.toStr(cnvFiles, ",");
        }
        if (plinkRoots != null) {
            processCommand += " plink=" + Array.toStr(plinkRoots, ",");
        }
        processCommand += " pop=" + popFile + " subPop=" + subPopFile;
        if (log1 != null) {
            processCommand += " log=" + log1.getFilename();
        }
        Files.write(processCommand, runDir + "../processResults.sh");
        Files.chmod(runDir + "../processResults.sh");
    }
    
    private static boolean checkRoots(String[] cnvFiles, String[] plinkRoots) {
        HashSet<String> rootSet = new HashSet<String>();
        if (plinkRoots != null) {
            for (String root : plinkRoots) {
                if (!rootSet.add(root)) {
                    return false;
                }
            }
        }
        if (cnvFiles != null) {
            for (String cnv : cnvFiles) {
                if (!rootSet.add(ext.rootOf(cnv, true))) {
                    return false;
                }
            }
        }
        return true;
    }

    private static void process(String dir, String[] cnvFiles, String[] plinkRoots, String popFile, String subPopFile, double thresh, Logger log1) throws IOException {
        PrintWriter writer;
        Logger log = log1 == null ? new Logger() : log1;
        String finalDir = ext.verifyDirFormat(dir);
        
        
        PopFileData popData, subPopData;
        try {
            popData = loadPopFile(popFile);
            subPopData = loadPopFile(subPopFile);
        } catch (IOException e) {
            e.printStackTrace();
            log.reportException(e);
            return;
        }
        
        String finalOut = finalDir + "final_results_pVals.xln";
        writer = Files.getAppropriateWriter(finalOut);
        boolean first = true;
        if (cnvFiles != null) {
            for (String cnvFile : cnvFiles) {
                String cnvRoot = ext.rootOf(cnvFile, true);
                String cnvDir = finalDir + cnvRoot + "/";
                writeResults(cnvDir + cnvRoot + "_results_pVals.xln", writer, cnvRoot, "#N/A", "#N/A", thresh, first);
                first = false;
                
                for (int p = 0; p < popData.pops.length; p++) {
                    String popDir = cnvDir + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "/";
                    writeResults(popDir + cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_results_pVals.xln", writer, cnvRoot, popData.pops[p], "#N/A", thresh, first);
                    
                    for (int sP = 0; sP < subPopData.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "/";
                        writeResults(subPopDir + cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_" + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "_results_pVals.xln", writer, cnvRoot, popData.pops[p], subPopData.pops[sP], thresh, first);
                    }
                }
            }
        }
        if (plinkRoots != null) {
            for (String plinkRoot : plinkRoots) {
                String plinkDir = finalDir + plinkRoot + "/";
                writeResults(plinkDir + plinkRoot + "_results_pVals.xln", writer, plinkRoot, "#N/A", "#N/A", thresh, first);
                first = false;
                
                for (int p = 0; p < popData.pops.length; p++) {
                    String popDir = plinkDir + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "/";
                    writeResults(popDir + plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_results_pVals.xln", writer, plinkRoot, popData.pops[p], "#N/A", thresh, first);
                    
                    for (int sP = 0; sP < subPopData.pops.length; sP++) {
                        String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "/";
                        writeResults(subPopDir + plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true) + "_" + ext.replaceWithLinuxSafeCharacters(subPopData.pops[sP], true) + "_results_pVals.xln", writer, plinkRoot, popData.pops[p], subPopData.pops[sP], thresh, first);
                    }
                }
            }
        }
        
        writer.flush();
        writer.close();
    }
    
    private static void writeResults(String file, PrintWriter writer, String root1, String root2, String root3, double thresh, boolean includeHeader) throws IOException {
        BufferedReader reader = Files.getAppropriateReader(file);
        
        String line = includeHeader ? null : reader.readLine();
        int pValueColumn = includeHeader ? -1 : ext.indexOfStr("tdt_P", line.split("\t"));
        boolean isHeader = true;
        while ((line = reader.readLine()) != null) {
            if (includeHeader && isHeader) {
                isHeader = false; 
                writer.println("DataFileRoot\tPopulation\tSubPopulation\t" + line);
                pValueColumn = ext.indexOfStr("tdt_P", line.split("\t"));
            } else {
                String val = line.split("\t")[pValueColumn];
                if (!ext.isMissingValue(val) && Double.parseDouble(val) < thresh) {
                    writer.println(root1 + "\t" + root2 + "\t" + root3 + "\t" + line);
                }
            }
        }
        reader.close();
    }

    private static void writeQsubs(ArrayList<String> qsubFiles, String runDir, String qsubQueue) {
        ArrayList<String> qsubCommands = new ArrayList<String>();
        for (String qsub : qsubFiles) {
            
            String dir = ext.parseDirectoryOfFile(qsub);
            qsubCommands.add("cd " + dir);
            qsubCommands.add("qsub " + (qsubQueue != null ? "-q " + qsubQueue : "") + ext.rootOf(qsub, true) + ".pbs");
            String ret = "cd ";
            for (int i = 0; i < dir.split("/").length; i++) {
                ret += "../";
            }
            qsubCommands.add(ret);
            
        }

        Files.writeArrayList(qsubCommands, runDir + "runPBSFiles.sh");
        Files.chmod(runDir + "runPBSFiles.sh");
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String runDir = "./";
        String[] cnvFiles = null;
        String[] plinkRoots = null;
        String pedFile = "./pedigree.dat";
        String popFile = "./pops.xln";
        String subPopFile = "./subPops.xln";
        double pThreshold = 0.05;
        String logFile = null;
        Logger log = null;
        String qsub = null;
        boolean process = false;
        
        String usage = "\n" + 
                       "gwas.EmimPipeline requires 2-4 arguments\n" + 
                       "   (1) run directory (i.e. dir=" + runDir + " (default))\n" + 
                       " AND\n" +
                       "   (2a) cnv files (i.e. cnvs=cnvFile1.cnv,cnvFile2.cnv (not the default))\n" +
                       " AND/OR \n" + 
                       "   (2b) PLINK fileroots (i.e. plink=plink1,plink2 (not the default))\n" + 
                       " AND, if desired (though the script to run this will be created automatically)\n" + 
                       "   (3) -process flag to consolidate results after PBS files have completed (i.e. -process (not the default))\n";
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("dir=")) {
                runDir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("cnvs=")) {
                cnvFiles = args[i].split("=")[1].split(",");
                numArgs--;
            } else if (args[i].startsWith("plink=")) {
                plinkRoots = args[i].split("=")[1].split(",");
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                logFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("pop=")) {
                popFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("subPop=")) {
                subPopFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("pThreshold=")) {
                pThreshold = ext.parseDoubleArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith("-process")) {
                process = true;
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (args.length == 0 || numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            if (logFile != null) {
                log = new Logger(logFile);
            }
            if (process) {
                process(runDir, cnvFiles, plinkRoots, popFile, subPopFile, pThreshold, log);
            } else {
                setup(runDir, cnvFiles, plinkRoots, pedFile, popFile, subPopFile, pThreshold, qsub, log);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
}
