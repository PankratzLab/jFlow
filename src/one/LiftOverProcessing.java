package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import common.CmdLine;
import common.Files;
import common.PSF;
import common.ext;

public class LiftOverProcessing {
    
    private static final String EXEC_LOC = "/home/pankrat2/public/bin/";
    private static final String EXEC_NAME = "liftOver";
    
    private static final String[] DB_LOCS = {
        "/home/pankrat2/public/bin/ref/hg18ToHg19.over.chain",
    };
    
    private static final boolean FORCE_BINARY = true;
    
    private static final String TEMP_EXT = ".posTemp";
    private static final String RESULT_EXT = ".outTemp";
    private static final String UNMAPPED_EXT = ".unTemp";
    private static final String DROPS_EXT = ".drops";
    private static final String TEMP_BIM = ".bimTemp";
    private static final String NEWFILES_APPEND = "_lifted";
    
    
    protected static void prepPlinkFiles(String plinkDirAndRoot, int conversion, boolean deleteTempFiles) {
        boolean writeTemp, writeDrops, writeNew;
        String cmd;
        HashMap<String, int[]> results;
        
        if (conversion >= DB_LOCS.length) {
            System.err.println("Error - invalid conversion");
            return;
        }
        
        if (!PSF.Plink.allFilesExist(plinkDirAndRoot, FORCE_BINARY)) {
            System.err.println("Error - unable to find all binary PLINK files (bed/bim/fam) in given directory with specified root: " + plinkDirAndRoot);
            return;
        }
        
        System.out.println("## Writing temporary positions file...");
        try {
            writeTemp = writeTempPositions(plinkDirAndRoot);
        } catch (IOException e) {
            e.printStackTrace();
            writeTemp = false;
        }
        if (!writeTemp) {
            return;
        }
        
        cmd = EXEC_LOC + EXEC_NAME + " " + 
                    plinkDirAndRoot + TEMP_EXT + " " + 
                    DB_LOCS[conversion] + " " + 
                    plinkDirAndRoot + RESULT_EXT + " " + 
                    plinkDirAndRoot + UNMAPPED_EXT;
        System.out.println("## Running liftOver...");
        CmdLine.runDefaults(cmd, "./");
        
        if (!Files.exists(plinkDirAndRoot + RESULT_EXT)) {
            System.err.println("Error - liftOver failed to output results.");
            return;
        }
        if (!Files.exists(plinkDirAndRoot + UNMAPPED_EXT)) {
            System.err.println("WARNING - no unmapped result file found - either liftOver failed to run, or all markers were successfully mapped.");
        }
        if (deleteTempFiles) {
            (new File(plinkDirAndRoot + TEMP_EXT)).delete();
        }
        
        System.out.println("## Parsing markers to remove due to lack of conversion...");
        try {
            writeDrops = parseDrops(plinkDirAndRoot);
        } catch (IOException e) {
            e.printStackTrace();
            writeDrops = false;
        }
        if (!writeDrops) {
            return;
        }
        if (deleteTempFiles) {
            (new File(plinkDirAndRoot + UNMAPPED_EXT)).delete();
        }
        
        System.out.println("## Excluding unconverted markers...");
        cmd = "plink2 --bfile " + plinkDirAndRoot + " --exclude " + plinkDirAndRoot + DROPS_EXT + " --make-bed --out " + plinkDirAndRoot + NEWFILES_APPEND;
        CmdLine.runDefaults(cmd, "./");
        if (deleteTempFiles) {
            (new File(plinkDirAndRoot + DROPS_EXT)).delete();
        }
        
        System.out.println("## Loading liftOver results...");
        try {
            results = loadResults(plinkDirAndRoot);
        } catch (IOException e) {
            e.printStackTrace();
            results = null;
        }
        if (results == null) {
            return;
        }
        if (deleteTempFiles) {
            (new File(plinkDirAndRoot + RESULT_EXT)).delete();
        }
        
        System.out.println("## Moving new .bim file to a temp file...");
        (new File(PSF.Plink.getBIM(plinkDirAndRoot + NEWFILES_APPEND))).renameTo(new File(plinkDirAndRoot + NEWFILES_APPEND + TEMP_BIM));
        
        System.out.println("## Writing new .bim file...");
        try {
            writeNew = writeNew(plinkDirAndRoot + NEWFILES_APPEND + TEMP_BIM, PSF.Plink.getBIM(plinkDirAndRoot + NEWFILES_APPEND), results);
        } catch (IOException e) {
            e.printStackTrace();
            writeNew = false;
        }
        if (!writeNew) {
            return;
        }
        if (deleteTempFiles) {
            (new File(plinkDirAndRoot + NEWFILES_APPEND + TEMP_BIM)).delete();
        }
        
        System.out.println("liftOver processing complete!");
    }
    
    private static boolean writeNew(String bimIn, String bimOut, HashMap<String, int[]> results) throws IOException {
        BufferedReader reader;
        PrintWriter writer, writerMissing;
        String line;
        String temp[];
        int missing;
        
        reader = Files.getAppropriateReader(bimIn);
        writer = Files.getAppropriateWriter(bimOut);
        writerMissing = Files.getAppropriateWriter(bimOut + ".missing");
        missing = 0;
        line = null;
        while ((line = reader.readLine()) != null) {
            temp = line.split("\t");
            int[] data = results.get(temp[1]);
            if (data == null) {
                missing++;
                writerMissing.println(line);
                writer.println(line);
                continue;
            }
            writer.println(data[0] + "\t" + temp[1] + "\t0\t" + data[1] + "\t" + temp[4] + "\t" + temp[5]);
        }
        writer.flush();
        writerMissing.flush();
        writer.close();
        writerMissing.close();
        reader.close();
        
        if (missing == 0) {
            (new File(bimOut + ".missing")).delete();
        }
        
        return true;
    }

    private static HashMap<String, int[]> loadResults(String plinkDirAndRoot) throws IOException {
        BufferedReader reader;
        HashMap<String, int[]> results;
        String line;
        String[] temp;
        
        reader = Files.getAppropriateReader(plinkDirAndRoot + RESULT_EXT);
        results = new HashMap<String, int[]>();
        line = null;
        while ((line = reader.readLine()) != null) {
            temp = line.split("\t");
            results.put(temp[3], new int[]{Integer.parseInt(temp[0].substring(3)), Integer.parseInt(temp[2])});
        }
        reader.close();
        
        return results;
    }
    
    private static boolean parseDrops(String plinkDirAndRoot) throws IOException {
        BufferedReader reader;
        PrintWriter writer;
        String line;
        String[] temp;
        
        reader = Files.getAppropriateReader(plinkDirAndRoot + UNMAPPED_EXT);
        writer = Files.getAppropriateWriter(plinkDirAndRoot + DROPS_EXT);
        
        line = null;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }
            temp = line.split("\t");
            writer.println(temp[3]);
        }
        writer.flush();
        writer.close();
        reader.close();
        return true;
    }
    
    private static boolean writeTempPositions(String plinkDirAndRoot) throws IOException {
        BufferedReader reader;
        PrintWriter writer;
        String tempFile, line;
        String[] temp;
        
        tempFile = plinkDirAndRoot + TEMP_EXT;
        if (Files.exists(tempFile)) {
            System.err.println("Error - temporary positions file already exists - did liftOver previously fail to complete?  Please remove and re-run.");
            return false;
        }
        
        reader = Files.getAppropriateReader(PSF.Plink.getBIM(plinkDirAndRoot));
        writer = Files.getAppropriateWriter(tempFile); 
        
        line = null;
        while ((line = reader.readLine()) != null) {
            temp = line.split("\t");
            writer.println("chr" + temp[0] + "\t" + (Integer.parseInt(temp[3]) - 1) + "\t" + temp[3] + "\t" + temp[1]);
        }
        writer.flush();
        writer.close();
        reader.close();
        
        return true;
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = "./plink";
        int conversion = 0; // hg18 to hg19
        boolean deleteTemp = true;

        String usage = "\n" + 
                       "one.ben.LiftOverProcessing requires 1+ arguments\n" + 
                       "   (1) PLINK directory and root (i.e. plink=" + filename + " (default))\n" + 
                       "   (2) delete temporary files? (i.e. deleteTemp=" + deleteTemp + " (default))\n" + 
                       "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("plink=")) {
                filename = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("deleteTemp=")) {
                deleteTemp = ext.parseBooleanArg(args[i]);
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
            prepPlinkFiles(filename, conversion, deleteTemp);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
