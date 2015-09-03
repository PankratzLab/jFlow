package one;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import common.Files;
import common.ext;

public class PLINK2GenomePackager {
    
//    plink2 --dosage fullList.txt list format=1 Zout --fam gedi_exome_plink.fam --covar GEDI_covars.dat --pheno GEDI_pheno_mtPC0_exome.dat
    
    void setup(String dir, String fileList) {
        
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
        
        ArrayList<String> plinkRuns = new ArrayList<String>();
        
        for (String pheno : phenos) {
            String[] parts = pheno.split("_");
            if (parts.length == 1) {
                System.err.println(ext.getTime() + "]\tError - phenotype file {" + pheno + "} is named incorrectly.  Correct naming is *_pheno.dat for covariate files.");
                continue;
            } else if (parts.length == 2) {
                if (covarPrefices.containsKey(parts[0])) {
                    plinkRuns.add(createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno));
                } else {
                    System.err.println(ext.getTime() + "]\tError - no covariate file found for phenotype file {" + pheno + "}.");
                    continue;
                }
            } else if (parts.length > 2) {
                String phenoCovarCheck = pheno.substring(0, pheno.length() - "_pheno.dat".length());
                if (covarPrefices.containsKey(phenoCovarCheck)) {
                    plinkRuns.add(createScript(fileList, famFile, covarPrefices.get(phenoCovarCheck), pheno));
                } else if (covarPrefices.containsKey(parts[0])) {
                    plinkRuns.add(createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno));
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
        
        
    }
    
    boolean setupForScript(String dir, String fileList, String famFile, String covarFile, String phenoFile) {
        
        String phenoPrefix = phenoFile.substring(0, phenoFile.length() - "_pheno.dat".length());
        
        String phenoDir = dir + phenoPrefix + "/";
        File phenoDirFile = new File(phenoDir);
        if (phenoDirFile.exists()) {
            return false;
        }
        boolean createdPhenoDir;
        try {
            createdPhenoDir = phenoDirFile.createNewFile();
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
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            return false;
        }
        
        return true;
    }
    
    String createScript(String fileList, String famFile, String covarFile, String phenoFile) {
        String script = "plink2 --dosage ../" + fileList + " list format=1 Zout --fam ../" + famFile + " --covar " + covarFile + " --pheno " + phenoFile;
        return script;
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String dir = "D:/test/";
        String fileList = "fullList.txt";

        String usage = "\n"
                        + "one.PLINK2GenomePackager requires 0-1 arguments\n"
                        + "   (1) directory (i.e. dir=" + dir + " (default))\n"
                        + "   (2) file list of source data files, indexed (i.e. file=" + fileList + " (default))\n"
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
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            new PLINK2GenomePackager().setup(dir, fileList);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
