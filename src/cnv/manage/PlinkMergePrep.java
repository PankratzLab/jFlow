package cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;

import common.Array;
import common.Files;
import common.HashVec;
import common.Positions;
import common.ext;

public class PlinkMergePrep {
    
    private static final int INVALID_ROOT = 0;
    private static final int PEDMAP = 1;
    private static final int BEDBIMFAM = 2;
    
    static int isValidRoot(String dir, String plinkRoot) {
        String[] files = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith(plinkRoot);
            }
        });
        boolean foundBed = false;
        boolean foundBim = false;
        boolean foundFam = false;
        boolean foundPed = false;
        boolean foundMap = false;
        for (String f : files) {
            if (f.endsWith(".bim")) {
                foundBim = true;
            } else if (f.endsWith(".bed")) {
                foundBed = true;
            } else if (f.endsWith(".fam")) {
                foundFam = true;
            } else if (f.endsWith(".ped")) {
                foundPed = true;
            } else if (f.endsWith(".map")) {
                foundMap = true;
            } 
        }
        
        if (foundBed && foundBim && foundFam) {
            return BEDBIMFAM;
        } else if (foundPed && foundMap) {
            return PEDMAP;
        }
        return INVALID_ROOT;
    }
    
    static String writeTempRegionFileFromUCSCFile(String regionFile, String dir) {
        String outFile = "tempRegions.txt";
        int index = 1;
        while (Files.exists(dir + outFile)) {
            outFile = "tempRegions_" + index++ + ".txt"; 
        }
        String[] ucscLines = HashVec.loadFileToStringArray(dir + regionFile, false, new int[]{0}, false);
        String[] plinkLines = new String[ucscLines.length];
        for (int i = 0; i < ucscLines.length; i++) {
            int[] pos = Positions.parseUCSClocation(ucscLines[i]);
            plinkLines[i] = pos[0] + "\t" + pos[1] + "\t" + pos[2] + "\tPLINK_" + i; 
        }
        Files.writeList(plinkLines, dir + outFile);
        return outFile;
    }
    
    static void renameMarkers(String dir, String plinkRoot, int rootType) {
        BufferedReader reader;
        PrintWriter writer;
        String line, ext;
        String[] parts;
        
        ext = rootType == PEDMAP ? ".map" : ".bim";
        boolean moved = (new File(dir + plinkRoot + ext)).renameTo(new File(dir + plinkRoot + "_orig" + ext));
        if (moved) {
            try {
                reader = Files.getAppropriateReader(dir + plinkRoot + "_orig" + ext);
                writer = Files.getAppropriateWriter(dir + plinkRoot + ext);
                while ((line = reader.readLine()) != null) {
                    parts = line.split("[\\s]+", -1);
                    parts[1] = plinkRoot + "_" + parts[1];
                    writer.println(Array.toStr(parts, "\t"));
                }
                writer.flush();
                writer.close();
                reader.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            System.err.println("Error - unable to move " + plinkRoot + ext + " to " + plinkRoot + "_orig" + ext + ".");
        }
    }
    
    static String merge(int outType, String outRoot, boolean overwrite, boolean renameMarkers, String regionFile, String dir, String plinkRoot1, String plinkRoot2) {
        int rootType1, rootType2;
        StringBuilder cmds;
        String tempRgnFile;
        
        if (outType != PEDMAP && outType != BEDBIMFAM) {
            System.err.println("Error - output type must be either (" + PEDMAP + ") PED/MAP, or (" + BEDBIMFAM + ") BED/BIM/FAM.");
            return null;
        }
        rootType1 = isValidRoot(dir, plinkRoot1);
        if (rootType1 == INVALID_ROOT) {
            System.err.println("Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '" + plinkRoot1 + "'.");
            return null;
        }
        rootType2 = isValidRoot(dir, plinkRoot2);
        if (rootType2 == INVALID_ROOT) {
            System.err.println("Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '" + plinkRoot2 + "'.");
            return null;
        }
        if (Files.exists(dir + outRoot + (outType == PEDMAP ? ".ped" : ".bed"))) {
            if (!overwrite) {
                System.err.println("Error - PLINK files with given output root ('" + outRoot + "') already exist!  Please specify alternate output name or set overwrite flag.");
                return null;
            }
        }
        if (regionFile != null && !regionFile.equals("") && !Files.exists(regionFile)) {
            System.err.println("Error - specified region file ('" + regionFile + "') doesn't exist!");
            return null;
        }
        dir = ext.verifyDirFormat(dir);
        
        if (renameMarkers) {
            renameMarkers(dir, plinkRoot1, rootType1);
            renameMarkers(dir, plinkRoot2, rootType2);
        }
        
        cmds = new StringBuilder();
        cmds.append("plink2 --noweb --").append(rootType1 == PEDMAP ? "file " : "bfile ").append(plinkRoot1).append(" --merge ");
        if (rootType2 == PEDMAP) {
            cmds.append(plinkRoot2).append(".ped ").append(plinkRoot2).append(".map ");
        } else if (rootType2 == BEDBIMFAM) {
            cmds.append(plinkRoot2).append(".bed ").append(plinkRoot2).append(".bim ").append(plinkRoot2).append(".fam ");
        }
        cmds.append("--").append(outType == PEDMAP ? "recode" : "make-bed");
        if (regionFile != null && !regionFile.equals("") && Files.exists(regionFile)) {
            tempRgnFile = writeTempRegionFileFromUCSCFile(regionFile, dir);
            cmds.append(" --extract range ").append(tempRgnFile);
        }
        cmds.append(" --out ").append(outRoot);
        return cmds.toString();
    }
    
    static String merge(int outType, String outRoot, boolean overwrite, boolean renameMarkers, String regionFile, String dir, String[] plinkRoots) {
        int[] rootTypes;
        String[] lines;
        String fileListName, tempRgnFile;
        StringBuilder cmds;
        int index;
        
        if (outType != PEDMAP && outType != BEDBIMFAM) {
            System.err.println("Error - output type must be either (" + PEDMAP + ") PED/MAP, or (" + BEDBIMFAM + ") BED/BIM/FAM.");
            return null;
        }
        if (plinkRoots.length < 2) {
            System.err.println("Error - at least two distinct PLINK file roots must be specified for merging.");
            return null;
        } else if (plinkRoots.length == 2) {
            return merge(outType, outRoot, overwrite, renameMarkers, regionFile, dir, plinkRoots[0], plinkRoots[1]);
        }
        if (Files.exists(dir + outRoot + (outType == PEDMAP ? ".ped" : ".bed"))) {
            if (!overwrite) {
                System.err.println("Error - PLINK files with given output root ('" + outRoot + "') already exist!  Please specify alternate output name or set overwrite flag.");
                return null;
            }
        }
        if (regionFile != null && !regionFile.equals("") && !Files.exists(dir + regionFile)) {
            System.err.println("Error - specified region file ('" + regionFile + "') doesn't exist!");
            return null;
        }
        dir = ext.verifyDirFormat(dir);
        
        rootTypes = new int[plinkRoots.length];
        for (int i = 0; i < plinkRoots.length; i++) {
            rootTypes[i] = isValidRoot(dir, plinkRoots[i]);
            if (rootTypes[i] == INVALID_ROOT) {
                System.err.println("Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '" + plinkRoots[i] + "'.");
                return null;
            }
        }
        lines = new String[plinkRoots.length - 1];
        for (int i = 1; i < plinkRoots.length; i++) {
            lines[i - 1] = rootTypes[i] == PEDMAP ? plinkRoots[i] + ".ped " + plinkRoots[i] + ".map" : plinkRoots[i] + ".bed " + plinkRoots[i] + ".bim " + plinkRoots[i] + ".fam";  
        }
        fileListName = "fileList.txt";
        index = 1;
        while (Files.exists(dir + fileListName)) {
            fileListName = "fileList_" + index++ + ".txt"; 
        }
        Files.writeList(lines, dir + fileListName);
        
        if (renameMarkers) {
            for (int i = 0; i < plinkRoots.length; i++) {
                renameMarkers(dir, plinkRoots[i], rootTypes[i]);
            }
        }
        
        cmds = new StringBuilder();
        cmds.append("plink2 --noweb --").append(rootTypes[0] == PEDMAP ? "file " : "bfile ").append(plinkRoots[0]).append(" --merge-list ").append(fileListName);
        cmds.append(" --").append(outType == PEDMAP ? "recode" : "make-bed");
        if (regionFile != null && !regionFile.equals("") && Files.exists(dir + regionFile)) {
            tempRgnFile = writeTempRegionFileFromUCSCFile(regionFile, dir);
            cmds.append(" --extract range ").append(tempRgnFile);
        }
        cmds.append(" --out ").append(outRoot);
        
        return cmds.toString();
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String outRoot = "merge";
        int outType = BEDBIMFAM;
        boolean overwrite = false;
        boolean renameMarkers = false;
        String dir = "./";
        String rgnFile = "";
        String[] roots = null;
        
        String usage = "\n" + 
                       "cnv.manage.PlinkMergePrep requires 0-1 arguments\n" + 
                       "   (1) Output Root (i.e. outRoot=" + outRoot + " (default))\n" +
                       "   (2) Output Type (i.e. outType=" + outType + " (default, options are (" + PEDMAP + ") PED/MAP, or (" + BEDBIMFAM + ") BED/BIM/FAM.))\n" +
                       "   (3) Region file (UCSC-format) (i.e. regions=" + rgnFile + " (default))\n" + 
                       "   (4) directory of files (i.e. dir=" + dir + " (default))\n" + 
                       "   (5) List of Plink fileRoots to merge (i.e. roots=plink1,plink2,plink4 (not the default))\n" + 
                       "   (6) overwrite flag (i.e. -overwrite (not the default))\n" + 
                       "   (7) rename markers flag (i.e. -renameMarkers (not the default))\n" + 
                       "";

        if (numArgs == 0) {
            System.err.println(usage);
            System.exit(1);
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("outRoot=")) {
                outRoot = ext.parseStringArg(args[i], "merge");
                numArgs--;
            } else if (args[i].startsWith("outType=")) {
                outType = ext.parseIntArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith("dir=")) {
                dir = ext.parseStringArg(args[i], "./");
                numArgs--;
            } else if (args[i].startsWith("roots=")) {
                roots = ext.parseStringArg(args[i], "").split(",");
                numArgs--;
            } else if (args[i].startsWith("regions=")) {
                rgnFile = ext.parseStringArg(args[i], "");
                numArgs--;
            } else if (args[i].startsWith("-overwrite")) {
                overwrite = true;
                numArgs--;
            } else if (args[i].startsWith("-renameMarkers")) {
                renameMarkers = true;
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
            String mergeCommand = merge(outType, outRoot, overwrite, renameMarkers, rgnFile, dir, roots);
            System.out.println(mergeCommand);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
