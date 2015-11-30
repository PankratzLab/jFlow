package cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import cnv.manage.SourceFileParser;
import common.Array;
import common.Elision;
import common.Files;
import common.Logger;
import common.ext;

public class SourceFileHeaderData implements Serializable {
    /**
     * 
     */
    private static final long serialVersionUID = -6906302109843776908L;
    
    String gsgtVersion = null;
    String processingDate = null;
    String content = null;
    int numSnps = -1;
    int totalSnps = -1;
    int numSamples = -1;
    int totalSamples = -1;
    int numFiles = -1;
    int currFile = -1;
    public int columnHeaderLineIndex = -1;
    
    public int colSampleIdent = -1;
    public int colSnpIdent = -1;
    public int colGenoAB1 = -1;
    public int colGenoAB2 = -1;
    public int colGeno1 = -1;
    public int colGeno2 = -1;
    public int colX = -1;
    public int colY = -1;
    public int colTheta = -1;
    public int colR = -1;
    public int colXRaw = -1;
    public int colYRaw = -1;
    public int colBAF = -1;
    public int colLRR = -1;
    public int colGC = -1;
    public String headerString = "";
    public String delimiter = "";
    public String[] cols = new String[0];
    
//    Order from ParseIllumina
//  0    GC
//  1    XRAW
//  2    YRAW
//  3    X
//  4    Y
//  5    Theta
//  6    R
//  7    B Allele Freq
//  8    Log R Ratio
    
    private SourceFileHeaderData() {}
    
    public static SourceFileHeaderData parseHeader(String file, Logger log) throws Elision, IOException {
        BufferedReader reader = Files.getAppropriateReader(file);
        String line = null;
        int lineCnt = 0;
        SourceFileHeaderData frhd = new SourceFileHeaderData();
        String delim = ",";
        while ((line = reader.readLine()) != null) {
            delim = ext.determineDelimiter(line); // TODO if file ends with .csv [or contains .csv?], can assume that delim is ','.  Sim., if ends with .xln, can assume delim is '\t'
            if (",".equals(delim)) {
                delim = "[\\s]*,[\\s]*"; // ext.indexFactors doesn't call trim()
            }
            if("[Data]".equals(line) || line.startsWith("rs") || line.toUpperCase().startsWith("SNP") || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true, false, false)[0] != -1) {
                break;
            }
            String[] parts = line.trim().split(",");
            processLine(parts, frhd);
            lineCnt++;
        }
        if (!"[Data]".equals(line) && !(line.startsWith("rs") || line.toUpperCase().startsWith("SNP") || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true, false, false)[0] != -1)) {
            log.reportError("Error - malformed or missing header.");
            throw new Elision(file);
        }
        if ("[Data]".equals(line)) {
            line = reader.readLine();
            delim = file.contains(".csv") ? "[\\s]*,[\\s]*" : file.contains(".xln") ? "[\\s]*\t[\\s]*" : ext.determineDelimiter(line);
            lineCnt++;
            if (!(line.startsWith("rs") || line.toUpperCase().startsWith("SNP") || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true, false, false)[0] != -1)) {
                log.reportError("Error - malformed or missing header.  Header must start with 'rs' or 'SNP' or contain one of the following: " + Array.toStr(SourceFileParser.SNP_HEADER_OPTIONS[0]) + ".");
                throw new Elision(file);
            }
        }
        reader.close(); // done
        reader = null; // release
        String columnHeaders = line;
        delim = file.contains(".csv") ? SOURCE_FILE_DELIMITERS.COMMA.getDelimiter() : file.contains(".xln") ? SOURCE_FILE_DELIMITERS.TAB.getDelimiter() : SOURCE_FILE_DELIMITERS.getDelimiter(ext.determineDelimiter(columnHeaders)).getDelimiter();
//        if (",".equals(delim)) {
//            delim = "[\\s]*,[\\s]*";
//        }
        parseColumnsBestGuess(columnHeaders.split(delim), frhd);
        frhd.setSourceFileDelimiter(delim);
        frhd.setHeaderString(columnHeaders);
        frhd.columnHeaderLineIndex = lineCnt;
        return frhd;
    }

    private static final String[] SAMPLE_FIELD_ID = { "Sample ID", "Sample Name" };
    
    private static final String[] DATA_FIELDS_GC = {"GC Score", "GCscore", "confidence", "confidenceScore"}; 
    private static final String[] DATA_FIELDS_XRAW = {"X Raw"};
    private static final String[] DATA_FIELDS_YRAW = {"Y Raw"}; 
    private static final String[] DATA_FIELDS_X = {"X", "Xvalue", "Log Ratio", "intensity_1", "Signal A"}; 
    private static final String[] DATA_FIELDS_Y = {"Y", "Yvalue", "Strength", "intensity_2", "Signal B"};
    private static final String[] DATA_FIELDS_THETA = {"Theta"};
    private static final String[] DATA_FIELDS_R = {"R"};
    private static final String[] DATA_FIELDS_BAF = {"B Allele Freq", "BAF"}; 
    private static final String[] DATA_FIELDS_LRR = {"Log R Ratio", "LRR"};
    private static final String[] GENOTYPE_FIELDS_A1_FOR = {"Allele1 - Forward", "Allele1", "genotype1", "Allele1 - Top", "Forward Strand Base Calls", "Forced Call", "Forced Call Codes"};
    private static final String[] GENOTYPE_FIELDS_A2_FOR = {"Allele2 - Forward", "Allele B", "genotype2", "Allele2 - Top", "Forward Strand Base Calls", "Forced Call", "Forced Call Codes"}; 
    private static final String[] GENOTYPE_FIELDS_A1_AB = {"Allele1 - AB", "Call Codes", "Call"}; 
    private static final String[] GENOTYPE_FIELDS_A2_AB = {"Allele2 - AB", "Call Codes", "Call"};
    
    private static final String[][] LOOKUP = {
/*0*/   SourceFileParser.SNP_HEADER_OPTIONS[0],
/*1*/   SAMPLE_FIELD_ID,
/*2*/   DATA_FIELDS_GC,
/*3*/   DATA_FIELDS_XRAW,
/*4*/   DATA_FIELDS_YRAW, 
/*5*/   DATA_FIELDS_X,
/*6*/   DATA_FIELDS_Y,
/*7*/   DATA_FIELDS_THETA,
/*8*/   DATA_FIELDS_R,
/*9*/   DATA_FIELDS_BAF,
/*10*/  DATA_FIELDS_LRR,
/*11*/  GENOTYPE_FIELDS_A1_FOR,
/*12*/  GENOTYPE_FIELDS_A2_FOR, 
/*13*/  GENOTYPE_FIELDS_A1_AB, 
/*14*/  GENOTYPE_FIELDS_A2_AB,
    };
    
    private static void parseColumnsBestGuess(String[] parts, SourceFileHeaderData frhd) throws Elision {
        int[] indices = ext.indexFactors(LOOKUP, parts, false, true, false, false);
        if (indices[0] == -1) {
            throw new Elision("Error - missing SNP ID column");
        }
        frhd.cols = parts;
        frhd.colSnpIdent = indices[0];
        frhd.colSampleIdent = indices[1];
//        frhd.col_sampleIndex = indices[2];
        frhd.colGC = indices[2];
        frhd.colXRaw = indices[3];
        frhd.colYRaw = indices[4];
        frhd.colX = indices[5];
        frhd.colY = indices[6];
        frhd.colTheta = indices[7];
        frhd.colR = indices[8];
        frhd.colBAF = indices[9];
        frhd.colLRR = indices[10];
        frhd.colGeno1 = indices[11];
        frhd.colGeno2 = indices[12];
        frhd.colGenoAB1 = indices[13];
        frhd.colGenoAB2 = indices[14];
    }
    
    private static void processLine(String[] parts, SourceFileHeaderData frhd) {
        if ("[Header]".equals(parts[0])) return;
        if ("File".equals(parts[0])) {
            String[] fileParts = parts[parts.length - 1].split(" of ");
            if (ext.isValidInteger(fileParts[0])) {
                frhd.currFile = Integer.parseInt(fileParts[0]);
            } else {
                // TODO error
            }
            if (ext.isValidInteger(fileParts[1])) {
                frhd.numFiles = Integer.parseInt(fileParts[1]);
            } else {
                // TODO error
            }
        }
        if ("Total Samples".equals(parts[0])) {
            if (ext.isValidInteger(parts[parts.length - 1])) {
                frhd.totalSamples = Integer.parseInt(parts[parts.length - 1]);
            } else {
                // TODO error
            }
        }
        if ("Num Samples".equals(parts[0])) {
            if (ext.isValidInteger(parts[parts.length - 1])) {
                frhd.numSamples = Integer.parseInt(parts[parts.length - 1]);
            } else {
                // TODO error
            }
        }
        if ("Total SNPs".equals(parts[0])) {
            if (ext.isValidInteger(parts[parts.length - 1])) {
                frhd.totalSnps = Integer.parseInt(parts[parts.length - 1]);
            } else {
                // TODO error
            }
        }
        if ("Num SNPs".equals(parts[0])) {
            if (ext.isValidInteger(parts[parts.length - 1])) {
                frhd.numSnps = Integer.parseInt(parts[parts.length - 1]);
            } else {
                // TODO error
            }
        }
        if ("Content".equals(parts[0])) {
            frhd.content = parts[parts.length - 1];
        }
        if ("Processing Date".equals(parts[0])) {
            frhd.processingDate = parts[parts.length - 1];
        }
        if ("GSGT Version".equals(parts[0])) {
            frhd.gsgtVersion = parts[parts.length - 1];
        }
    }
    
    public static HashMap<String, SourceFileHeaderData> validate(final String rawDir, final String ext, boolean fullValidation, Logger log) {
        String dir = rawDir.endsWith("/") || rawDir.endsWith("\\") ? rawDir : common.ext.verifyDirFormat(rawDir);
        String[] possibleFiles = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(ext);
            }
        });
        
        boolean valid = false;
        HashMap<String, SourceFileHeaderData> headers = null;
        try {
//            if (fullValidation) {
                headers = new HashMap<String, SourceFileHeaderData>();
//            }
            for (String possFile : possibleFiles) {
                SourceFileHeaderData frhd = SourceFileHeaderData.parseHeader(dir + possFile, log);
//                if (fullValidation) {
                    headers.put(possFile, frhd);
//                }
            }
            if (fullValidation) {
                String error = doFullValidation(headers, log);
                if (error != null) {
                    throw new Elision("Error - " + error);
                }
            }
            valid = true;
        } catch (Elision e) {
            log.reportException(e);
        } catch (IOException e) {
            log.reportException(e);
        }
        
        return valid ? headers : null;
    }
    
    public static String doFullValidation(HashMap<String, SourceFileHeaderData> headers, Logger log) {
        int cnt = headers.size();
        HashMap<Integer, ArrayList<String>> numSnpsSet = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> totSnpsSet = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> numSamplesSet = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> totSamplesSet = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> headerLineIndex = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> sampleID = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> snpIndex = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> genoAB_1 = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> genoAB_2 = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> genoForward_1 = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> genoForward_2 = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> x = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> y = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> theta = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> r = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> xRaw = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> yRaw = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> BAF = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> LRR = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> gc = new HashMap<Integer, ArrayList<String>>();
        for (java.util.Map.Entry<String, SourceFileHeaderData> entry : headers.entrySet()) {
            SourceFileHeaderData headerData = entry.getValue();
            if (headerData.numFiles == -1) {
                headerData.numFiles = cnt;
            } else if (headerData.numFiles != cnt) {
                return "Number of Files listed in Source File {" + entry.getKey() + "} does not equal the number of headers needing validation.  Please check source directory and extension and try again.";
            }
            ArrayList<String> files = numSnpsSet.get(headerData.numSnps);
            if (files == null) {
                files = new ArrayList<String>();
                numSnpsSet.put(headerData.numSnps, files);
            }
            files.add(entry.getKey());
            files = totSnpsSet.get(headerData.totalSnps);
            if (files == null) {
                files = new ArrayList<String>();
                totSnpsSet.put(headerData.totalSnps, files);
            }
            files.add(entry.getKey());
            files = numSamplesSet.get(headerData.numSamples);
            if (files == null) {
                files = new ArrayList<String>();
                numSamplesSet.put(headerData.numSamples, files);
            }
            files.add(entry.getKey());
            files = totSamplesSet.get(headerData.totalSamples);
            if (files == null) {
                files = new ArrayList<String>();
                totSamplesSet.put(headerData.totalSamples, files);
            }
            files.add(entry.getKey());
            files = headerLineIndex.get(headerData.columnHeaderLineIndex);
            if (files == null) {
                files = new ArrayList<String>();
                headerLineIndex.put(headerData.columnHeaderLineIndex, files);
            }
            files.add(entry.getKey());
            files = sampleID.get(headerData.colSampleIdent);
            if (files == null) {
                files = new ArrayList<String>();
                sampleID.put(headerData.colSampleIdent, files);
            }
            files.add(entry.getKey());
//            files = sampleIndex.get(entry.getValue().col_sampleIndex);
//            if (files == null) {
//                files = new ArrayList<String>();
//                sampleIndex.put(entry.getValue().col_sampleIndex, files);
//            }
            files.add(entry.getKey());
            files = snpIndex.get(headerData.colSnpIdent);
            if (files == null) {
                files = new ArrayList<String>();
                snpIndex.put(headerData.colSnpIdent, files);
            }
            files.add(entry.getKey());
            files = genoAB_1.get(headerData.colGenoAB1);
            if (files == null) {
                files = new ArrayList<String>();
                genoAB_1.put(headerData.colGenoAB1, files);
            }
            files.add(entry.getKey());
            files = genoAB_2.get(headerData.colGenoAB2);
            if (files == null) {
                files = new ArrayList<String>();
                genoAB_2.put(headerData.colGenoAB2, files);
            }
            files.add(entry.getKey());
            files = genoForward_1.get(headerData.colGeno1);
            if (files == null) {
                files = new ArrayList<String>();
                genoForward_1.put(headerData.colGeno1, files);
            }
            files.add(entry.getKey());
            files = genoForward_2.get(headerData.colGeno2);
            if (files == null) {
                files = new ArrayList<String>();
                genoForward_2.put(headerData.colGeno2, files);
            }
            files.add(entry.getKey());
            files = x.get(headerData.colX);
            if (files == null) {
                files = new ArrayList<String>();
                x.put(headerData.colX, files);
            }
            files.add(entry.getKey());
            files = y.get(headerData.colY);
            if (files == null) {
                files = new ArrayList<String>();
                y.put(headerData.colY, files);
            }
            files.add(entry.getKey());
            files = theta.get(headerData.colTheta);
            if (files == null) {
                files = new ArrayList<String>();
                theta.put(headerData.colTheta, files);
            }
            files.add(entry.getKey());
            files = r.get(headerData.colR);
            if (files == null) {
                files = new ArrayList<String>();
                r.put(headerData.colR, files);
            }
            files.add(entry.getKey());
            files = xRaw.get(headerData.colXRaw);
            if (files == null) {
                files = new ArrayList<String>();
                xRaw.put(headerData.colXRaw, files);
            }
            files.add(entry.getKey());
            files = yRaw.get(headerData.colYRaw);
            if (files == null) {
                files = new ArrayList<String>();
                yRaw.put(headerData.colYRaw, files);
            }
            files.add(entry.getKey());
            files = BAF.get(headerData.colBAF);
            if (files == null) {
                files = new ArrayList<String>();
                BAF.put(headerData.colBAF, files);
            }
            files.add(entry.getKey());
            files = LRR.get(headerData.colLRR);
            if (files == null) {
                files = new ArrayList<String>();
                LRR.put(headerData.colLRR, files);
            }
            files.add(entry.getKey());
            files = gc.get(headerData.colGC);
            if (files == null) {
                files = new ArrayList<String>();
                gc.put(headerData.colGC, files);
            }
            files.add(entry.getKey());
        }
        int numErrors = 0;
        HashSet<String> errorFiles = new HashSet<String>();
        if (numSnpsSet.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : numSnpsSet.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (totSnpsSet.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : totSnpsSet.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (numSamplesSet.size() > 1) {
            numErrors++;            
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : numSamplesSet.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (totSamplesSet.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : totSamplesSet.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (headerLineIndex.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : headerLineIndex.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (sampleID.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : sampleID.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
//        if (sampleIndex.size() > 1) {
//            numErrors++;
//            ArrayList<String> minSet = null;
//            for (ArrayList<String> set : sampleIndex.values()) {
//                if (minSet == null || set.size() < minSet.size()) {
//                    minSet = set;
//                }
//            }
//            errorFiles.addAll(minSet);
//        }
        if (snpIndex.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : snpIndex.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (genoAB_1.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : genoAB_1.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (genoAB_2.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : genoAB_2.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (genoForward_1.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : genoForward_1.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (genoForward_2.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : genoForward_2.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (x.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : x.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (y.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : y.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (theta.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : theta.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (r.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : r.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (xRaw.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : xRaw.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (yRaw.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : yRaw.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (BAF.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : BAF.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (LRR.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : LRR.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        if (gc.size() > 1) {
            numErrors++;
            ArrayList<String> minSet = null;
            for (ArrayList<String> set : gc.values()) {
                if (minSet == null || set.size() < minSet.size()) {
                    minSet = set;
                }
            }
            errorFiles.addAll(minSet);
        }
        
        if (numErrors > 0 || errorFiles.size() > 0) {
            if (errorFiles.size() > 0) {
                for (String file : errorFiles) {
                    log.reportError("Error - data or column mismatch suspected in source file {" + file + "}.  Please verify file integrity and try again.");
                }
            }
            return "Found " + numErrors + " data or column index mismatches among " + errorFiles.size() + " files.  Please check log for more details";
        }
        
        return null;
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String dir = "D:/data/ny_registry/gedi_gwas/00src/";//null;
        String ext = ".csv.gz";
        String log = null;

        String usage = "\n" + 
                       "cnv.filesys.FinalReportHeaderData requires 2 argument\n" + 
                       "   (1) Directory of FinalReport files (i.e. dir=" + dir + " (default))\n" + 
                       "   (2) Extension of FinalReport files (i.e. ext=" + ext + " (default))\n";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("dir=")) {
                dir = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("ext=")) {
                ext = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                log = args[i].split("=")[1];
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0 || dir == null) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            validate(dir, ext, true, log == null ? new Logger() : new Logger(log));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String getHeaderString() {
        return headerString;
    }

    public void setHeaderString(String headerString) {
        this.headerString = headerString;
    }
    
    private void setSourceFileDelimiter(String delim) {
        this.delimiter = delim;
    }
    
    public String getSourceFileDelimiter() {
        return this.delimiter;
    }
}