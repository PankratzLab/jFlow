package one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.Vector;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Matrix;
import common.ext;

public class DBGAPMerge {
    
    Logger log;
    
    public static class FileSet {
        String dir;
        String dataFile;
        String dataDict;
        String varReport;
        
        ArrayList<DataColumn> dataDefs;
        ArrayList<String> ids = new ArrayList<String>();
        HashMap<String, String[]> idDataMap = new HashMap<String, String[]>();
    }
    
    public static class DataColumn {
        String source;
        String study;
        String table;
        String consentGroup;
        String varID;
        String varName;
        String varDesc;
        HashMap<String, String> varValueDescriptions = new HashMap<String, String>();
        String varUnit;
        String comment;
    }
    
    private static final String DICT_SUFF = "data_dict.xml";
    private static final String VAR_RPT_SUFF = "var_report.xml";
    private static final String DATA_SUFF = "txt.gz";
    
    private static final String MISSING_DATA = ".";
    
    public ArrayList<FileSet> discoverFiles(String dir) {
        ArrayList<FileSet> fss = new ArrayList<DBGAPMerge.FileSet>();
        
        String[] allFiles = new File(dir).list();
        
        for (String gzFile : allFiles) {
            if (!gzFile.endsWith(DATA_SUFF)) continue;
            String[] parts = gzFile.split("\\.");
            
            String key = parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];
            
            String dict = null;
            String varRpt = null;
            for (String f : allFiles) {
                if (f.endsWith(".gz") || !f.startsWith(key)) continue;
                if (f.endsWith(DICT_SUFF)) {
                    if (dict == null) {
                        dict = f;
                    } else {
                        System.err.println("Error - duplicate data dictionary file detected: " + dict + " / " + f);
                    }
                } else if (f.endsWith(VAR_RPT_SUFF)) {
                    if (varRpt == null) {
                        varRpt = f;
                    } else {
                        System.err.println("Error - duplicate variable report file detected: " + varRpt + " / " + f);
                    }
                }
            }
            
            if (dict != null && varRpt != null) {
                FileSet fs = new FileSet();
                fs.dir = dir;
                fs.dataFile = gzFile;
                fs.dataDict = dict;
                fs.varReport = varRpt;
                fss.add(fs);
            } else {
                if (dict == null) {
                    log.reportError("Error - no data dictionary file found for data file: " + gzFile);
                }
                if (varRpt == null) {
                    log.reportError("Error - no variable report file found for data file: " + gzFile);
                }
            }
        }
        
        return fss;
    }
    
    public ArrayList<DataColumn> readDataDict(FileSet fs) throws ParserConfigurationException, SAXException, IOException {
        ArrayList<DataColumn> dict = new ArrayList<DBGAPMerge.DataColumn>();
        
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse(new File(fs.dir + fs.dataDict));
        doc.getDocumentElement().normalize();
        
        Node tableNode = doc.getElementsByTagName("data_table").item(0);
        String table = tableNode.getAttributes().getNamedItem("id").getNodeValue();
        String study = tableNode.getAttributes().getNamedItem("study_id").getNodeValue();
        
        NodeList variables = doc.getElementsByTagName("variable");
        
        for (int i = 0, count = variables.getLength(); i < count; i++) {
            Node varNode = variables.item(i);
            
            DataColumn dc = new DataColumn();
            Node vnName = varNode.getAttributes().getNamedItem("id");
            dc.source = fs.dir + fs.dataFile;
            dc.varID = vnName.getNodeValue().trim();
            dc.table = table.trim();
            dc.study = study.trim();
            
            String[] pts = fs.dataFile.split("\\.");
            for (String s : pts) {
                if (s.length() == 2 && s.startsWith("c")) {
                    dc.consentGroup = s.trim().replaceAll("\\.", "");
                    break;
                }
            }
            
            NodeList nlChildren = varNode.getChildNodes();
            for (int c = 0, countCh = nlChildren.getLength(); c < countCh; c++) {
                Node chNode = nlChildren.item(c);
                String chNm = chNode.getNodeName();
                if (chNm.equals("#text")) {
                    continue;
                }
                String nodeNm = chNode.getNodeName();
                String nodeVal = chNode.getTextContent();
                if (nodeNm.equals("name")) {
                    dc.varName = nodeVal;
                } else if (nodeNm.equals("description")) {
                    dc.varDesc = nodeVal;
                } else if (nodeNm.equals("type")) {
                    
                } else if (nodeNm.equals("value")) {
                    NamedNodeMap nnm = chNode.getAttributes();
                    Node cdNode = nnm.getNamedItem("code");
                    if (cdNode != null) {
                        dc.varValueDescriptions.put(cdNode.getNodeValue(), nodeVal);
                    } else {
                        System.out.println("Value without Code: " + nodeVal);
                    }
                } else if (nodeNm.equals("unit")) {
                    dc.varUnit = nodeVal;
                } else if (nodeNm.equals("comment")) {
                    dc.comment = nodeVal.replaceAll("\\n", " - ");
                } else {
                    log.reportError("Unknown node: " + nodeNm + " = " + nodeVal);
                }
                
            }
            dict.add(dc);
        }
        
        return dict;
    }
    
    public void readDataFile(FileSet fs) throws IOException {
        String line;
        String[] parts;
        BufferedReader reader = Files.getAppropriateReader(fs.dir + fs.dataFile);
        
        boolean firstAfterHeader = false;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("##")) {
                // TODO read header line
                firstAfterHeader = true;
            } else if (line.startsWith("#") || line.equals("")) {
                continue; // skip unimportant header content
            } else {
                if (firstAfterHeader) {
                    firstAfterHeader = false;
                    continue;
                }
                parts = line.split("\t", -1);
                fs.ids.add(parts[0]);
                fs.idDataMap.put(parts[0], Array.subArray(parts, 1));
            }
        }
        reader.close();
        reader = null;
        System.gc();
    }
    
    public void run(String[] dirs, String outFile, String outMap, Logger log) {
        this.log = log;
        ArrayList<FileSet> files = new ArrayList<DBGAPMerge.FileSet>();
        for (String dir : dirs) {
            files.addAll(discoverFiles(ext.verifyDirFormat(dir)));
        }
        
        for (FileSet fs : files) {
            try {
                fs.dataDefs = readDataDict(fs);
            } catch (ParserConfigurationException e) {
                log.reportException(e);
            } catch (SAXException e) {
                log.reportException(e);
            } catch (IOException e) {
                log.reportException(e);
            }
            if (fs.dataDefs != null) {
                try {
                    readDataFile(fs);
                } catch (IOException e) {
                    log.reportException(e);
                }
            }
        }
        
        TreeSet<String> idSet = new TreeSet<String>();
        for (FileSet fs : files) {
            idSet.addAll(fs.ids);
        }
        
        PrintWriter writer = Files.getAppropriateWriter(outFile);
        StringBuilder lineOut = null;

        lineOut = new StringBuilder("DBGAP_ID");
        for (FileSet fs : files) {
            for (DataColumn dc : fs.dataDefs) {
                lineOut.append("\t").append(dc.varName).append(";").append(dc.varID).append(";").append(dc.table);
                if (dc.consentGroup != null) {
                    lineOut.append(";").append(dc.consentGroup);
                }
            }
        }
        writer.println(lineOut.toString());
        
        for (String id : idSet) {
            lineOut = new StringBuilder(id);
            for (FileSet fs : files) {
                String[] data = fs.idDataMap.get(id);
                if (data == null) {
                    data = Array.stringArray(fs.dataDefs.size(), MISSING_DATA);
                }
                for (String s : data) {
                    lineOut.append("\t").append("".equals(s) ? MISSING_DATA : s);
                }
            }
            writer.println(lineOut.toString());
        }
        writer.flush();
        writer.close();
        
        writer = Files.getAppropriateWriter(outMap);
        writer.println("Source\tStudy\tTable\tVariable\tVariableID\tfinalColumnName\tcustomColumnName\tDescription\tUnits\tVariableMapping\tComment");
        for (FileSet fs : files) {
            for (DataColumn dc : fs.dataDefs) {
                lineOut = new StringBuilder();
                lineOut.append(dc.source).append("\t");
                lineOut.append(dc.study).append("\t");
                lineOut.append(dc.table).append("\t");
                lineOut.append(dc.varName).append("\t");
                lineOut.append(dc.varID).append("\t");
                lineOut.append(dc.varName).append(";").append(dc.varID).append(";").append(dc.table);
                if (dc.consentGroup != null) {
                    lineOut.append(";").append(dc.consentGroup);
                }
                lineOut.append("\t");
                lineOut.append(".").append("\t");
                lineOut.append(dc.varDesc).append("\t");
                lineOut.append(dc.varUnit == null ? "." : dc.varUnit).append("\t");
                if (dc.varValueDescriptions.isEmpty()) {
                    lineOut.append(".").append("\t");
                } else {
                    for (Entry<String, String> ent : dc.varValueDescriptions.entrySet()) {
                        lineOut.append(ent.getKey()).append("=").append(ent.getValue()).append(";");
                    }
                }
                if (dc.comment != null && !"".equals(dc.comment)) {
                    lineOut.append("\t").append(dc.comment);
                } else {
                    lineOut.append("\t").append(".");
                }
                writer.println(lineOut.toString());
            }
        }
        writer.flush();
        writer.close();
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String out = null;
        String outMap = null;
        String[] dir = null;
        String logfile = null;
        Logger log;
        
        String usage = "\n" + 
                        "one.ben.DBGAPMergeAndLookup requires 3+ arguments\n" +
                        "To MERGE dbGap files:\n" + 
                        "   (1) Data Output Filename (i.e. out=" + out + " (default))\n" + 
                        "   (2) Map Output Filename (i.e. outMap=" + outMap + " (default))\n" + 
                        "   (3) Input directory (or a comma-delimited list of directories) (i.e. dir=" + dir + " (default))\n" + 
                        "   (4) OPTIONAL: Log file (i.e. log=" + logfile + " (default))\n" + 
                        "\n" + 
                        "";

//        boolean test = true;
//        if (test) {
//            DBGapLookup.fromParameters("F:/dbGap merge/search.crf", new Logger());
            
//            DBGapExtract.fromParameters("F:/dbGap merge/dbgap.crf", new Logger());
            
//            String fileDir1 = "/home/pankrat2/shared/ARIC/phenos/dbGaP/50859/PhenoGenotypeFiles/RootStudyConsentSet_phs000280.ARIC_RootStudy.v3.p1.c1.HMB-IRB/PhenotypeFiles/";
//            String fileDir2 = "/home/pankrat2/shared/ARIC/phenos/dbGaP/50865/PhenoGenotypeFiles/RootStudyConsentSet_phs000280.ARIC_RootStudy.v3.p1.c2.DS-CVD-IRB/PhenotypeFiles/";
//            dir = new String[]{fileDir1, fileDir2};
//            out = "/scratch.global/coleb/mergeOut.xln.gz";
//            outMap = "/scratch.global/coleb/mergeMap.xln";
//
//            log = new Logger(logfile);
//            (new DBGAPMerge()).run(dir, out, outMap, log);
//            return;
//        }
        
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("out=")) {
                out = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("outMap=")) {
                outMap = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("dir=")) {
                dir = ext.parseStringArrayArg(args[i], null);
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                logfile = args[i].split("=")[1];
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
            log = new Logger(logfile);
            (new DBGAPMerge()).run(dir, out, outMap, log);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static class DBGapExtract {
        
        public static void fromParameters(String filename, Logger log) {
            Vector<String> params;
            String[] args;
            
            params = Files.parseControlFile(filename, "dbgap", new String[]{
                                "variables=searchTerms.xln",
                                "dataFile=mergeOut.xln.gz",
                                "outputFile=searchedVariables.xln",
                                "head=IID",
                        }, log);
            
            if (params != null) {
                args = params.toArray(new String[params.size()]);
                main(args);
            }
        }
        
        public static void main(String[] args) {
            int numArgs = args.length;
            String varFile = "searchTerms.xln";
            String dataFile = "mergeOut.xln.gz";
            String outputFile = "searchedVariables.xln";
            String headIdent = "IID";
            String logfile = null;
            Logger log;

            String usage = "\n" + 
                            "one.ben.DBGAPExtract requires 4+ arguments\n" + 
                            "   (1) Searched variables file (output from DBGapLookup) (i.e. variables=" + varFile + " (default))\n" + 
                            "   (2) Merged dbGap data file (data output from DBGapMerge ) (i.e. dataFile=" + dataFile + " (default))\n" + 
                            "   (3) Extracted output filename (i.e. outputFile=" + outputFile + " (default))\n" + 
                            "   (4) ID column name in output file (i.e. head=" + headIdent + " (default))\n" +
                            "   (5) OPTIONAL: Log file (i.e. log=" + logfile + " (default))\n" + 
                            "";

            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                    System.err.println(usage);
                    System.exit(1);
                } else if (args[i].startsWith("variables=")) {
                    varFile = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("dataFile=")) {
                    dataFile = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("outputFile=")) {
                    outputFile = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("head=")) {
                    headIdent = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("log=")) {
                    logfile = args[i].split("=")[1];
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
                log = new Logger(logfile);
                extract(varFile, dataFile, outputFile, headIdent, log);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        private static void extract(String varFile, String dataFile, String outputFile, String headIdent, Logger log) throws IOException {
            BufferedReader dataReader;
            PrintWriter writer;
            String[][] varData;
            String[] colsToLoad, parts;
            String line, delim;
            StringBuilder sb;
            
            varData = HashVec.loadFileToStringMatrix(varFile, true, new int[]{1, 2}, false); // ident col, repl col, ignore source and other cols
            colsToLoad = Matrix.extractColumn(varData, 0);
            
            dataReader = Files.getAppropriateReader(dataFile);
            line = dataReader.readLine();
            delim = ext.determineDelimiter(line);
            parts = line.split(delim, -1);
            int[] factors = ext.indexFactors(colsToLoad, parts, false, false);
            
            writer = Files.getAppropriateWriter(outputFile);
            sb = new StringBuilder();
            sb.append(headIdent);
            for (int i = 0; i < factors.length; i++) {
                if (factors[i] >= 0) {
                    sb.append("\t").append(ext.isMissingValue(varData[i][1]) ? varData[i][0] : varData[i][1]);
                }
            }
            writer.println(sb.toString());
            
            while((line = dataReader.readLine()) != null) {
                parts = line.split(delim, -1);
                sb = new StringBuilder();
                sb.append(parts[0]);
                for (int i = 0; i < factors.length; i++) {
                    if (factors[i] >= 0) {
                        sb.append("\t").append(parts[factors[i]]);
                    }
                }
                writer.println(sb.toString());
            }
            writer.flush();
            writer.close();
            dataReader.close();
        }
        
    }
    
    public static class DBGapLookup {
        
        public static void fromParameters(String filename, Logger log) {
            Vector<String> params;
            String[] line, args;
            String mapFile;
            StringBuilder sb;
            
            params = Files.parseControlFile(filename, "search", new String[]{
                                                "mergeMap.xln   searchCols=3,7,9,10 outputCols=0,5,6,3,7,9,10 out=searchTerms.xln  crf=dbgap.crf",
                                                "# search terms below, one per line:"
                                            }, log);
            if (params != null) {
                line = params.remove(0).trim().split("[\\s]+");
                
                mapFile = line[0];
                if (!Files.exists(mapFile) && Files.exists(ext.verifyDirFormat(ext.parseDirectoryOfFile(filename)) + mapFile)) {
                    mapFile = ext.verifyDirFormat(ext.parseDirectoryOfFile(filename)) + mapFile;
                }

                args = new String[line.length + (log.getFilename() == null ? 1 : 2)];
                args[0] = "map=" + mapFile;
                for (int i = 1; i < line.length; i++) {
                    args[i] = line[i];
                }
                sb = new StringBuilder(params.get(0));
                for (int i = 1; i < params.size(); i++) {
                    sb.append(",").append(params.get(i));
                }
                args[line.length] = "search=" + sb.toString();
                if (log.getFilename() != null) {
                    args[line.length + 1] = "log=" + log.getFilename();
                }
                main(args);
            }
        }
        
        public static void main(String[] args) {
            int numArgs = args.length;
            
            String map = "mergeMap.xln";
            int[] searchCols = new int[]{3,7,9,10,11};
            int[] outputCols = new int[]{0,5,6,3,7,9,10,11};
            String output = "searchTerms.xln";
            String crf = "dbgap.crf";
            String[] search = new String[]{"cancer","height","bmi"};
            
            String logfile = null;
            Logger log;

            String usage = "\n" + 
                            "one.ben.DBGAPMerge requires 6+ arguments\n" + 
                            "   (1) dbGap Map File (from merged files) (i.e. map=" + map + " (default))\n" +
                            "   (2) Column indices in which to search (i.e. searchCols=" + Array.toStr(searchCols, ",") + " (default))\n" + 
                            "   (3) Column indices from which to export data (i.e. outputCols=" + Array.toStr(outputCols, ",") + " (default))\n" + 
                            "   (4) Output File (i.e. out=" + output + " (default))\n" +
                            "   (5) Data export CRF filename (i.e. crf=" + crf + " (default))\n" +
                            "   (6) Search terms (i.e. search=" + Array.toStr(search, ",") + " (default))\n" + 
                            "   (7) OPTIONAL: Log file (i.e. log=" + logfile + " (default))\n" + 
                            "";
            
            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                    System.err.println(usage);
                    System.exit(1);
                } else if (args[i].startsWith("map=")) {
                    map = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("crf=")) {
                    crf = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("out=")) {
                    output = args[i].split("=")[1];
                    numArgs--;
                } else if (args[i].startsWith("searchCols=")) {
                    searchCols = ext.parseIntArrayArg(args[i]);
                    numArgs--;
                } else if (args[i].startsWith("outputCols=")) {
                    outputCols = ext.parseIntArrayArg(args[i]);
                    numArgs--;
                } else if (args[i].startsWith("search=")) {
                    search = ext.parseStringArrayArg(args[i], search);
                    numArgs--;
                } else if (args[i].startsWith("log=")) {
                    logfile = args[i].split("=")[1];
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
                log = new Logger(logfile);
                search(map, search, searchCols, outputCols, output, crf, log);
                if (Files.exists(crf)) {
                    log.reportError("Error - specified .CRF file already exists! -- " + crf);
                } else {
                    Files.parseControlFile(crf, "dbgap", new String[]{
                            "variables=" + output,
                            "dataFile=mergeOut.xln.gz",
                            "outputFile=searchedVariables.xln",
                            "head=IID",
                    }, log);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        private static void search(String map, String[] search, int[] searchCols, int[] outputCols, String output, String crf, Logger log) throws IOException {
            String line, delim, outputDelim;
            String[] parts;
            StringBuilder sb;
            
            BufferedReader reader = Files.getAppropriateReader(map);
            PrintWriter writer = Files.getAppropriateWriter(output);
            
            line = reader.readLine().trim();
            delim = ext.determineDelimiter(line);
            outputDelim = "\t";
            parts = line.split(delim, -1);
            
            sb = new StringBuilder();
            for (int i = 0; i < outputCols.length; i++) {
                sb.append(parts[outputCols[i]]).append(outputDelim);
            }
            for (int i = 0; i < search.length; i++) {
                sb.append(search[i]);
                if (i < search.length - 1) {
                    sb.append(outputDelim);
                }
            }
            writer.println(sb.toString());
            
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if ("".equals(line)) continue;
                parts = line.trim().split(delim, -1);

                int[] searchInds = Array.intArray(search.length, -1);
                for (int s = 0; s < search.length; s++) {
                    for (int col : searchCols) {
                        searchInds[s] = Math.max(searchInds[s], parts[col].toLowerCase().indexOf(search[s].toLowerCase()));
                    }
                }
                boolean foundAny = false;
                for (int s : searchInds) {
                    if (s >= 0) {
                        foundAny = true;
                        break;
                    }
                }
                if (foundAny) {
                    sb = new StringBuilder();
                    for (int i = 0; i < outputCols.length; i++) {
                        sb.append(parts[outputCols[i]]).append(outputDelim);
                    }
                    for (int i = 0; i < searchInds.length; i++) {
                        sb.append(searchInds[i]);
                        if (i < search.length - 1) {
                            sb.append(outputDelim);
                        }
                    }
                    writer.println(sb.toString());
                }
            }
            reader.close();
            writer.flush();
            writer.close();
        }
        
    }
    
    
    
}
