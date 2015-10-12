package cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import common.Elision;
import common.Files;
import common.Logger;
import common.ext;

public class FinalReportHeaderData {
    String gsgtVersion = null;
    String processingDate = null;
    String content = null;
    int numSnps = -1;
    int totalSnps = -1;
    int numSamples = -1;
    int totalSamples = -1;
    int numFiles = 1;
    int currFile = 1;
    int headerLineIndex = -1;
    
    public int col_sampleID = -1;
//    int col_sampleIndex = -1;
//    int col_snpName = -1;
    public int col_snpIndex = -1;
    public int col_genoAB_1 = -1;
    public int col_genoAB_2 = -1;
    public int col_geno_1 = -1;
    public int col_geno_2 = -1;
    public int col_x = -1;
    public int col_y = -1;
    public int col_theta = -1;
    public int col_r = -1;
    public int col_xRaw = -1;
    public int col_yRaw = -1;
    public int col_BAF = -1;
    public int col_LRR = -1;
    public int col_gc = -1;
    private String headerString = "";
    
    public FinalReportHeaderData() {}
    
    public static FinalReportHeaderData parseHeader(String file) throws Elision, IOException {
        BufferedReader reader = Files.getAppropriateReader(file);
//        if (!"[Header]".equals(reader.readLine())) {
//            // TODO error, malformed header
//            throw new Elision(file);
//        }
        String line = null;
        int lineCnt = 1;
        FinalReportHeaderData frhd = new FinalReportHeaderData();
        while ((line = reader.readLine()) != null && !"[Data]".equals(line) && !(line.startsWith("rs") || line.toUpperCase().startsWith("SNP"))) {
            String[] parts = line.trim().split(",");
            processLine(parts, frhd);
            lineCnt++;
        }
        if (!"[Data]".equals(line)) {
            // TODO error, malformed header
            throw new Elision(file);
        }
        String columnHeaders = reader.readLine();
        reader.close(); // done
        parseColumnsBestGuess(columnHeaders.split("[\\s]*,[\\s]*"), frhd); // TODO double check regex for consuming whitespace padding
        frhd.setHeaderString(columnHeaders);
        frhd.headerLineIndex = lineCnt + 1;
        return frhd;
    }

    private static final String[] SNP_HEADER_OPTIONS = { "SNP Name", "rsID", "Probe Set ID", "ProbeSet" };
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
/*0*/   SNP_HEADER_OPTIONS,
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
    
    private static void parseColumnsBestGuess(String[] parts, FinalReportHeaderData frhd) throws Elision {
        int[] indices = ext.indexFactors(LOOKUP, parts, false, true, false, false);
        if (indices[0] == -1) {
            throw new Elision("Error - missing SNP ID column");
        }
        frhd.col_snpIndex = indices[0];
        frhd.col_sampleID = indices[1];
//        frhd.col_sampleIndex = indices[2];
        frhd.col_gc = indices[2];
        frhd.col_xRaw = indices[3];
        frhd.col_yRaw = indices[4];
        frhd.col_x = indices[5];
        frhd.col_y = indices[6];
        frhd.col_theta = indices[7];
        frhd.col_r = indices[8];
        frhd.col_BAF = indices[9];
        frhd.col_LRR = indices[10];
        frhd.col_geno_1 = indices[11];
        frhd.col_geno_2 = indices[12];
        frhd.col_genoAB_1 = indices[13];
        frhd.col_genoAB_2 = indices[14];
    }
    
    private static void processLine(String[] parts, FinalReportHeaderData frhd) {
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
    
    public static HashMap<String, FinalReportHeaderData> validate(final String dir, final String ext, boolean fullValidation, Logger log) {
        String[] possibleFiles = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(ext);
            }
        });
        
        boolean valid = false;
        HashMap<String, FinalReportHeaderData> headers = null;
        try {
            if (fullValidation) {
                headers = new HashMap<String, FinalReportHeaderData>();
            }
            for (String possFile : possibleFiles) {
                FinalReportHeaderData frhd = FinalReportHeaderData.parseHeader(dir + possFile);
                if (fullValidation) {
                    headers.put(possFile, frhd);
                }
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
    
    public static String doFullValidation(HashMap<String, FinalReportHeaderData> headers, Logger log) {
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
        for (java.util.Map.Entry<String, FinalReportHeaderData> entry : headers.entrySet()) {
            if (entry.getValue().numFiles != cnt) {
                return "Number of Files listed in Source File {" + entry.getKey() + "} does not equal the number of headers needing validation.  Please check source directory and extension and try again.";
            }
            ArrayList<String> files = numSnpsSet.get(entry.getValue().numSnps);
            if (files == null) {
                files = new ArrayList<String>();
                numSnpsSet.put(entry.getValue().numSnps, files);
            }
            files.add(entry.getKey());
            files = totSnpsSet.get(entry.getValue().totalSnps);
            if (files == null) {
                files = new ArrayList<String>();
                totSnpsSet.put(entry.getValue().totalSnps, files);
            }
            files.add(entry.getKey());
            files = numSamplesSet.get(entry.getValue().numSamples);
            if (files == null) {
                files = new ArrayList<String>();
                numSamplesSet.put(entry.getValue().numSamples, files);
            }
            files.add(entry.getKey());
            files = totSamplesSet.get(entry.getValue().totalSamples);
            if (files == null) {
                files = new ArrayList<String>();
                totSamplesSet.put(entry.getValue().totalSamples, files);
            }
            files.add(entry.getKey());
            files = headerLineIndex.get(entry.getValue().headerLineIndex);
            if (files == null) {
                files = new ArrayList<String>();
                headerLineIndex.put(entry.getValue().headerLineIndex, files);
            }
            files.add(entry.getKey());
            files = sampleID.get(entry.getValue().col_sampleID);
            if (files == null) {
                files = new ArrayList<String>();
                sampleID.put(entry.getValue().col_sampleID, files);
            }
            files.add(entry.getKey());
//            files = sampleIndex.get(entry.getValue().col_sampleIndex);
//            if (files == null) {
//                files = new ArrayList<String>();
//                sampleIndex.put(entry.getValue().col_sampleIndex, files);
//            }
            files.add(entry.getKey());
            files = snpIndex.get(entry.getValue().col_snpIndex);
            if (files == null) {
                files = new ArrayList<String>();
                snpIndex.put(entry.getValue().col_snpIndex, files);
            }
            files.add(entry.getKey());
            files = genoAB_1.get(entry.getValue().col_genoAB_1);
            if (files == null) {
                files = new ArrayList<String>();
                genoAB_1.put(entry.getValue().col_genoAB_1, files);
            }
            files.add(entry.getKey());
            files = genoAB_2.get(entry.getValue().col_genoAB_2);
            if (files == null) {
                files = new ArrayList<String>();
                genoAB_2.put(entry.getValue().col_genoAB_2, files);
            }
            files.add(entry.getKey());
            files = genoForward_1.get(entry.getValue().col_geno_1);
            if (files == null) {
                files = new ArrayList<String>();
                genoForward_1.put(entry.getValue().col_geno_1, files);
            }
            files.add(entry.getKey());
            files = genoForward_2.get(entry.getValue().col_geno_2);
            if (files == null) {
                files = new ArrayList<String>();
                genoForward_2.put(entry.getValue().col_geno_2, files);
            }
            files.add(entry.getKey());
            files = x.get(entry.getValue().col_x);
            if (files == null) {
                files = new ArrayList<String>();
                x.put(entry.getValue().col_x, files);
            }
            files.add(entry.getKey());
            files = y.get(entry.getValue().col_y);
            if (files == null) {
                files = new ArrayList<String>();
                y.put(entry.getValue().col_y, files);
            }
            files.add(entry.getKey());
            files = theta.get(entry.getValue().col_theta);
            if (files == null) {
                files = new ArrayList<String>();
                theta.put(entry.getValue().col_theta, files);
            }
            files.add(entry.getKey());
            files = r.get(entry.getValue().col_r);
            if (files == null) {
                files = new ArrayList<String>();
                r.put(entry.getValue().col_r, files);
            }
            files.add(entry.getKey());
            files = xRaw.get(entry.getValue().col_xRaw);
            if (files == null) {
                files = new ArrayList<String>();
                xRaw.put(entry.getValue().col_xRaw, files);
            }
            files.add(entry.getKey());
            files = yRaw.get(entry.getValue().col_yRaw);
            if (files == null) {
                files = new ArrayList<String>();
                yRaw.put(entry.getValue().col_yRaw, files);
            }
            files.add(entry.getKey());
            files = BAF.get(entry.getValue().col_BAF);
            if (files == null) {
                files = new ArrayList<String>();
                BAF.put(entry.getValue().col_BAF, files);
            }
            files.add(entry.getKey());
            files = LRR.get(entry.getValue().col_LRR);
            if (files == null) {
                files = new ArrayList<String>();
                LRR.put(entry.getValue().col_LRR, files);
            }
            files.add(entry.getKey());
            files = gc.get(entry.getValue().col_gc);
            if (files == null) {
                files = new ArrayList<String>();
                gc.put(entry.getValue().col_gc, files);
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
}