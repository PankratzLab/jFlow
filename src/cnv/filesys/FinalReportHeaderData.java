package cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;

import common.Elision;
import common.Files;
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
    
    int col_sampleID = -1;
    int col_sampleIndex = -1;
//    int col_snpName = -1;
    int col_snpIndex = -1;
    int col_genoAB_1 = -1;
    int col_genoAB_2 = -1;
    int col_genoForward_1 = -1;
    int col_genoForward_2 = -1;
    int col_x = -1;
    int col_y = -1;
    int col_theta = -1;
    int col_r = -1;
    int col_xRaw = -1;
    int col_yRaw = -1;
    int col_BAF = -1;
    int col_LRR = -1;
    int col_gc = -1;
    String headerString;
    
    private FinalReportHeaderData() {}
    
    public static FinalReportHeaderData parseHeader(String file) throws Elision, IOException {
        BufferedReader reader = Files.getAppropriateReader(file);
        if (!"[Header]".equals(reader.readLine())) {
            // TODO error, malformed header
            throw new Elision(file);
        }
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
        frhd.headerString = columnHeaders;
        frhd.headerLineIndex = lineCnt + 1;
        return frhd;
    }

    private static final String[] SNP_HEADER_OPTIONS = { "SNP Name", "rsID", "Probe Set ID", "ProbeSet" };
    private static final String[] SAMPLE_FIELD_ID = { "Sample ID", "Sample Name" };
    private static final String[] SAMPLE_FIELD_INDEX = { "Sample Index" };
    
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
/*2*/   SAMPLE_FIELD_INDEX,
/*3*/   DATA_FIELDS_GC,
/*4*/   DATA_FIELDS_XRAW,
/*5*/   DATA_FIELDS_YRAW, 
/*6*/   DATA_FIELDS_X,
/*7*/   DATA_FIELDS_Y,
/*8*/   DATA_FIELDS_THETA,
/*9*/   DATA_FIELDS_R,
/*10*/   DATA_FIELDS_BAF,
/*11*/  DATA_FIELDS_LRR,
/*12*/  GENOTYPE_FIELDS_A1_FOR,
/*13*/  GENOTYPE_FIELDS_A2_FOR, 
/*14*/  GENOTYPE_FIELDS_A1_AB, 
/*15*/  GENOTYPE_FIELDS_A2_AB,
    };
    
    private static void parseColumnsBestGuess(String[] parts, FinalReportHeaderData frhd) throws Elision {
        int[] indices = ext.indexFactors(LOOKUP, parts, false, true, false, false);
        if (indices[0] == -1) {
            throw new Elision("Error - missing SNP ID column");
        }
        frhd.col_snpIndex = indices[0];
        frhd.col_sampleID = indices[1];
        frhd.col_sampleIndex = indices[2];
        frhd.col_gc = indices[3];
        frhd.col_xRaw = indices[4];
        frhd.col_yRaw = indices[5];
        frhd.col_x = indices[6];
        frhd.col_y = indices[7];
        frhd.col_theta = indices[8];
        frhd.col_r = indices[9];
        frhd.col_BAF = indices[10];
        frhd.col_LRR = indices[11];
        frhd.col_genoForward_1 = indices[12];
        frhd.col_genoForward_2 = indices[13];
        frhd.col_genoAB_1 = indices[14];
        frhd.col_genoAB_2 = indices[15];
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
    
    public static FinalReportHeaderData validate(final String dir, final String ext, boolean fullValidation) {
        String[] possibleFiles = (new File(dir)).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(ext);
            }
        });
        
        boolean valid = false;
        FinalReportHeaderData returnValue = null;
        try {
            HashMap<String, FinalReportHeaderData> headers = null;
            if (fullValidation) {
                headers = new HashMap<String, FinalReportHeaderData>();
            }
            for (String possFile : possibleFiles) {
                FinalReportHeaderData frhd = FinalReportHeaderData.parseHeader(dir + possFile);
                if (returnValue == null) {
                    returnValue = frhd;
                }
                if (fullValidation) {
                    headers.put(possFile, frhd);
                }
            }
            if (fullValidation) {
                doFullValidation(headers);
            }
            valid = true;
        } catch (Elision e) {
            // TODO error, final report file malformed!
            e.printStackTrace();
        } catch (IOException e) {
            // TODO error reading file report file
            e.printStackTrace();
        }
        
        return valid ? returnValue : null;
    }
    
    private static void doFullValidation(HashMap<String, FinalReportHeaderData> headers) {
        int cnt = headers.size();
        for (java.util.Map.Entry<String, FinalReportHeaderData> entry : headers.entrySet()) {
            
        }
        // TODO verify all headers are the same - #s are all equal, etc
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String dir = "D:/data/ny_registry/gedi_gwas/00src/";//null;
        String ext = ".csv.gz";

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
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0 || dir == null) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            validate(dir, ext, true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}