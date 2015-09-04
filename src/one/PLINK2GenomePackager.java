package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import common.Files;
import common.HashVec;
import common.ext;

public class PLINK2GenomePackager {
    
//    plink2 --dosage fullList.txt list format=1 Zout --fam gedi_exome_plink.fam --covar GEDI_covars.dat --pheno GEDI_pheno_mtPC0_exome.dat
    
    void setup(String dir, String fileList, String pmAllFile) {
        
        File dirFile = new File(dir);
        
        if (!dirFile.exists()) {
            System.err.println("Error - run directory { " + dir + " } does not exist");
            return;
        }
        
        String[] covars = dirFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith("_covars.dat");
            }
        });
        if (covars.length == 0) {
            System.err.println(ext.getTime() + "]\tError - no covariate files found.  Ensure that files are named *_covars.dat and try again.");
            return;
        }
        
        
        HashMap<String, String> covarPrefices = new HashMap<String, String>();
        
        for (String covarFile : covars) {
            String[] parts = covarFile.split("_");
            if (parts.length == 1) {
                System.err.println(ext.getTime() + "]\tError - covariate file {" + covarFile + "} is named incorrectly.  Correct naming is *_covars.dat for covariate files.");
                continue;
            } else if (parts.length == 2) {
                covarPrefices.put(parts[0], covarFile);
            } else if (parts.length > 2) {
                covarPrefices.put(covarFile.substring(0, covarFile.length() - "_covars.dat".length()), covarFile);
            }
        }
        
        String[] phenos = dirFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith("_pheno.dat");
            }
        });
        if (phenos.length == 0) {
            System.err.println(ext.getTime() + "]\tError - no phenotype files found.  Ensure that files are named *_pheno.dat and try again.");
            return;
        }
        
        String[] famList = dirFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".fam");
            }
        });
        if (famList.length == 0) {
            System.err.println(ext.getTime() + "]\tError - no .fam file found.");
            return;
        } else if (famList.length > 1) {
            System.err.println(ext.getTime() + "]\tError - multiple .fam files found.");
            return;
        }
        
        String famFile = famList[0];
        
        HashMap<String, String> plinkRuns = new HashMap<String, String>();
        
        for (String pheno : phenos) {
            String[] parts = pheno.split("_");
            if (parts.length == 1) {
                System.err.println(ext.getTime() + "]\tError - phenotype file {" + pheno + "} is named incorrectly.  Correct naming is *_pheno.dat for covariate files.");
                continue;
            } else if (parts.length == 2) {
                if (covarPrefices.containsKey(parts[0])) {
                    if (setupForScript(dir, fileList, famFile, covarPrefices.get(parts[0]), pheno)) {
                        String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
                        String phenoDir = dir + phenoPrefix + "/";
                        String script = createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno, phenoDir); 
                        plinkRuns.put(phenoDir, script);
                    }
                } else {
                    System.err.println(ext.getTime() + "]\tError - no covariate file found for phenotype file {" + pheno + "}.");
                    continue;
                }
            } else if (parts.length > 2) {
                String phenoCovarCheck = pheno.substring(0, pheno.length() - "_pheno.dat".length());
                if (covarPrefices.containsKey(phenoCovarCheck)) {
                    if (setupForScript(dir, fileList, famFile, covarPrefices.get(phenoCovarCheck), pheno)) {
                        String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
                        String phenoDir = dir + phenoPrefix + "/";
                        String script = createScript(fileList, famFile, covarPrefices.get(phenoCovarCheck), pheno, phenoDir);
                        plinkRuns.put(phenoDir, script);
                    }
                } else if (covarPrefices.containsKey(parts[0])) {
                    if (setupForScript(dir, fileList, famFile, covarPrefices.get(parts[0]), pheno)) {
                        String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
                        String phenoDir = dir + phenoPrefix + "/";
                        String script = createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno, phenoDir);
                        plinkRuns.put(phenoDir, script);
                    }
                } else {
                    System.err.println(ext.getTime() + "]\tError - no covariate file found for phenotype file {" + pheno + "}.");
                    continue;
                }
            }
        }
        
        if (plinkRuns.size() == 0) {
            System.err.println(ext.getTime() + "]\tError - no valid PLINK2 scripts found.");
            return;
        }
        
        for (java.util.Map.Entry<String, String> plinkRun : plinkRuns.entrySet()) {
            Files.qsub(plinkRun.getKey() + "runPlink2.qsub", plinkRun.getValue(), 5000, 4, 1);
        }
        
        StringBuilder sb = new StringBuilder("cd ").append(dir).append("\n");
        for (String plinkDir : plinkRuns.keySet()) {
            sb.append("cd ").append(plinkDir).append("\n")
                .append("qsub runPlink2.qsub").append("\n")
                .append("cd ..").append("\n");
        }
        Files.write(sb.toString(), dir + "masterRun.sh");
        Files.chmod(dir + "masterRun.sh");
        
        setupProcess(dir, plinkRuns, pmAllFile);
    }
    
    boolean setupForScript(String dir, String fileList, String famFile, String covarFile, String phenoFile) {
        
        String phenoPrefix = phenoFile.substring(0, phenoFile.length() - "_pheno.dat".length());
        
        String phenoDir = dir + phenoPrefix + "/";
        File phenoDirFile = new File(phenoDir);
        if (phenoDirFile.exists()) {
            return false;
        }
        boolean createdPhenoDir;
        createdPhenoDir = phenoDirFile.mkdir();
        if (!createdPhenoDir) {
            System.err.println(ext.getTime() + "]\tError - couldn't create directory {" + phenoDir + "}.");
            return false;
        }
        boolean fileCopy = Files.copyFile(dir + covarFile, phenoDir + covarFile);
        if (!fileCopy) {
            System.err.println(ext.getTime() + "]\tError - couldn't copy covariate file into phenotype subdirectory {" + phenoDir + "}.");
            return false;
        }
        fileCopy = Files.copyFile(dir + phenoFile, phenoDir + phenoFile);
        if (!fileCopy) {
            System.err.println(ext.getTime() + "]\tError - couldn't copy phenotype file into phenotype subdirectory {" + phenoDir + "}.");
            return false;
        }
        fileCopy = Files.copyFile(dir + famFile, phenoDir + famFile);
        if (!fileCopy) {
            System.err.println(ext.getTime() + "]\tError - couldn't copy .fam file into phenotype subdirectory {" + phenoDir + "}.");
            return false;
        }
        
        return true;
    }
    
    String createScript(String fileList, String famFile, String covarFile, String phenoFile, String phenoDir) {
        String script = "cd " + phenoDir + "\nplink2 --dosage ../" + fileList + " list format=1 Zout --fam " + famFile + " --covar " + covarFile + " --pheno " + phenoFile;
        return script;
    }
    
    void setupProcess(String dir, HashMap<String, String> plinkRuns, String pmAllFile) {
        StringBuilder masterProc = new StringBuilder();
        for (String plinkDir : plinkRuns.keySet()) {
            StringBuilder procCmd = new StringBuilder();
            procCmd.append("cd " + plinkDir + "\n")
            .append("jcp one.PLINK2GenomePackager -process dir=")
            .append(plinkDir)
            .append(" pmFile=")
            .append(pmAllFile)
            .append(" N=`grep -o -E '[0-9]+ people pass filters and QC' *.o | sed 's/[^0-9]*//g'`")
            .append("\n");
            Files.qsub(plinkDir + "process.qsub", procCmd.toString(), 5000, 3, 1);
            masterProc.append("cd ").append(plinkDir).append("\n");
            masterProc.append("qsub process.qsub").append("\n");
        }
        Files.write(masterProc.toString(), dir + "masterProcess.sh");
        Files.chmod(dir + "masterProcess.sh");
    }
    
    void process(String dir, String pmAllFile, String N) {
        if (N == null || "".equals(N) || !ext.isValidInteger(N)) {
            System.err.println(ext.getTime() + "]\tError - specified value of N {" + N + "} is invalid.");
            return; // continue without N?
        }
        
        if (!(new File(pmAllFile)).exists()) {
            System.err.println(ext.getTime() + "]\tError - specified PM_ALL file {" + pmAllFile + "} doesn't exist.");
            return;
        }
        
        String file = dir + "plink.assoc.dosage.gz";
        String newHeader = "SNP\tCHR\tPOS\tA1\tA2\tN\tFRQ\tINFO\tBETA\tSE\tP";
        
        Hashtable<String, String> pmAll = HashVec.loadFileToHashString(pmAllFile, 1, new int[]{0, 3}, "\t", false, false);
        
        try {
            BufferedReader reader = Files.getAppropriateReader(file);
            PrintWriter writer = Files.getAppropriateWriter(dir + ext.rootOf(dir) + ".plink.assoc.dosage.gz");
            writer.println(newHeader);
            String line = null;
            reader.readLine();
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                String[] parts = line.split("[\\s]+");
                String mkrChrPos = pmAll.get(parts[0]);
                
                StringBuilder sb = new StringBuilder();
                sb.append(parts[0]).append("\t")
                    .append(parts[1]).append("\t")
                    .append(parts[2]).append("\t");
                if (mkrChrPos != null) {
                    sb.append(mkrChrPos);
                } else {
                    sb.append(".\t.");
                }
                for (int i = 3; i < parts.length; i++) {
                    sb.append("\t").append(parts[i]);
                }
                writer.println(sb.toString());
            }
            writer.flush();
            writer.close();
            reader.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
    }
    
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String dir = "D:/test/";
        String fileList = "fullList.txt";
        boolean process = false;
        String N = null;
        String pmAll = "D:/test/PM_all";
        
        String usage = "\n"
                        + "one.PLINK2GenomePackager requires 3-4 arguments\n"
                        + "   (1) directory (i.e. dir=" + dir + " (default))\n"
                        + "   (2) file list of source data files, indexed (i.e. file=" + fileList + " (default))\n"
                        + "   (3) location of PM_ALL file (i.e. pmFile=PM_ALL (not the default))\n"
                        + " OR:\n"
                        + "   (1) directory of PLINK2 run (i.e. dir=" + dir + " (not the default))\n"
                        + "   (2) location of PM_ALL file (i.e. pmFile=PM_ALL (not the default))\n"
                        + "   (3) Number of individuals in analysis (i.e. N=1000 (not the default))\n"
                        + "   (4) -process flag"
                        + "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("dir=")) {
                dir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("file=")) {
                fileList = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("pmFile=")) {
                pmAll = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("N=")) {
                N = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("-process")) {
                process = true;
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            if (process) {
                new PLINK2GenomePackager().process(dir, pmAll, N);
            } else {
                new PLINK2GenomePackager().setup(dir, fileList, pmAll);
            }   
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
