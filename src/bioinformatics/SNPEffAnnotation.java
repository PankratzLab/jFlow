package bioinformatics;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import common.Aliases;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEffCmdEff;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

public class SNPEffAnnotation {
    
    private static final String[][] FACTORS = new String[][]{Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS};
    private static final String[][] FACTORS_5 = new String[][]{Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS, Aliases.ALLELES[0], Aliases.ALLELES[1]};
    private static final int DEFAULT_BUILD = 37;
    private static final String DEFAULT_BUILD_STR = "hg19";
    
    public static String processInput(String infile, Logger log) {
        File f = new File(infile);
        if (!f.exists()) {
            log.reportError("Error - file {" + infile + "} doesn't exist.");
            return null; 
        }
        String[] firstLine = Files.getHeaderOfFile(infile, null);
        if (firstLine.length == 1) {
            // snps only
            if (ext.indexOfStr(firstLine[0], Aliases.MARKER_NAMES) == -1) {
                log.reportError("Error - files containing only one column of data must contain RS ID's");
                return null;
            }
            return prepareInput(infile, null, log);
            
        } else if (firstLine.length == 3) {
            // snps, chr, and pos
            int[] factors = ext.indexFactors(FACTORS, firstLine, false, true, true, false);
            if (Array.countIf(factors, -1) > 0) {
                log.reportError("Error - files containing three columns of data must contain RS ID's, Chromosomes, and Positions.");
                return null;
            }
            if (factors[0] == 2 && factors[1] == 0 && factors[2] == 1) {
                return infile;
            } else {
                return prepareInput(infile, factors, log);
            }
        } else if (firstLine.length >= 5) {
            // snps, chr, and pos
            int[] factors = ext.indexFactors(FACTORS_5, firstLine, false, true, true, false);
            if (Array.countIf(factors, -1) > 0) {
                log.reportError("Error - files containing five or more columns of data must contain RS ID's, Chromosomes, Positions, and Ref and Alt alleles.");
                return null;
            }
            if (factors[0] == 2 && factors[1] == 0 && factors[2] == 1) {
                return infile;
            } else {
                return prepareInput(infile, factors, log);
            }
        } else {
            log.reportError("Error - file {" + infile + "} must have 1 (rsID), 3 (rsID, chr, pos), or 5+ (rsID, chr, pos, ref, alt) columns.");
            return null;
        }
    }
    
    private static String prepareInput(String file, int[] factorIndices, Logger log) {
        String fileToUse = file;
        int[] indicesToUse = factorIndices;
        if (factorIndices == null) {
//            ParseSNPlocations.parseSNPlocations(snpListFile, vcfFile, unmappedVCF, mergedVCF, log, monitor);
            ParseSNPlocations.lowMemParse(file, MapSNPsAndGenes.getSNPDB(DEFAULT_BUILD, log), MapSNPsAndGenes.getMergeDB(log), true, log); // TODO using hash parse, not VCF
            fileToUse = ext.rootOf(file, false)+"_positions.xln";
            indicesToUse = ext.indexFactors(FACTORS, Files.getHeaderOfFile(fileToUse, log), false, true, true, false);
        }
        String newFile = ext.rootOf(file, false)+"_snpEffLookup.txt";
        int cnt = 0;
        while (new File(newFile).exists()) {
            newFile = ext.rootOf(file, false)+"_snpEffLookup_" + cnt++ + ".txt";
        }
        int[] colsToLoad = indicesToUse.length == 3 ? new int[]{indicesToUse[1], indicesToUse[2], indicesToUse[0]} : new int[]{indicesToUse[1], indicesToUse[2], indicesToUse[0], indicesToUse[3], indicesToUse[4]};
        String[] rschrpos = HashVec.loadFileToStringArray(fileToUse, false, true, colsToLoad, false);
        Files.writeList(rschrpos, newFile);
        return newFile;
    }
    
    public static void pipeline(String inputFile, String snpEffConfigFile, Logger log) {
        String snpEffInputFile = processInput(inputFile, log);
        String[] args = {"ann", "-v", "-c", snpEffConfigFile, DEFAULT_BUILD_STR, snpEffInputFile};
        // use -t [multithreading, implies -noStats]
        SnpEff snpEff = new SnpEff(args);
        SnpEffCmdEff snpEffCmd = (SnpEffCmdEff) snpEff.snpEffCmd(); // instance of SnpEffCmdEff
        List<VcfEntry> vcfList = snpEffCmd.run(true);
        PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(inputFile, false) + "_snpEff.out.vcf");
        for (VcfEntry vc : vcfList) {
            writer.println(vc.toString());
        }
        writer.flush();
        writer.close();
    }

    public static String getDefaultConfigFile() {
        return "./snpEff.config"; // TODO put snpEff.config file into valid jar location [test current location]
    }
    

    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = null;
        String config = SNPEffAnnotation.getDefaultConfigFile();
        String logFile = null;
        
        String usage = "\\n"+
        "bioinformatics.SNPEffAnnotation requires 1-3 arguments\n"+
        "   (1) name of file containing 1 (rsIDs), 3 (rsIDs, pos, chr), or 5 (rsIDs, pos, chr, ref, alt) data columns (i.e. file="+filename+" (default))\n"+
        "   (2) SNPEFF config file (i.e. config="+config+" (default))\n"+
        "   (3) log file (i.e. log=null (default))\n"+
        "";

        for (int i = 0; i<args.length; i++) {
            if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("file=")) {
                filename = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("config=")) {
                config = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                logFile = args[i].split("=")[1];
                numArgs--;
            }
        }
        if (numArgs != 0 || filename == null) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            SNPEffAnnotation.pipeline(filename, config, logFile == null ? new Logger() : new Logger(logFile));
        } catch (Exception e) {
            e.printStackTrace();      
        }
    }
    
}
