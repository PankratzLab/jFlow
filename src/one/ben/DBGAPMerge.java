package one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeSet;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import common.Array;
import common.Files;
import common.ext;

public class DBGAPMerge {
    
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
                    System.err.println("Error - no data dictionary file found for data file: " + gzFile);
                }
                if (varRpt == null) {
                    System.err.println("Error - no variable report file found for data file: " + gzFile);
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
        
        NodeList variables = doc.getElementsByTagName("variable");
        
        for (int i = 0, count = variables.getLength(); i < count; i++) {
            Node varNode = variables.item(i);
            
            DataColumn dc = new DataColumn();
            Node vnName = varNode.getAttributes().getNamedItem("id");
            dc.varID = vnName.getNodeValue();
            
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
                    dc.comment = nodeVal;
                } else {
                    
                    System.out.println("Unknown node: " + nodeNm + " = " + nodeVal);
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
    
    public void run(String[] dirs, String outFile, String outMap) {
        ArrayList<FileSet> files = new ArrayList<DBGAPMerge.FileSet>();
        for (String dir : dirs) {
            files.addAll(discoverFiles(ext.verifyDirFormat(dir)));
        }
        
        for (FileSet fs : files) {
            try {
                fs.dataDefs = readDataDict(fs);
            } catch (ParserConfigurationException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            } catch (SAXException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            if (fs.dataDefs != null) {
                try {
                    readDataFile(fs);
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
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
                lineOut.append("\t").append(dc.varName).append(";").append(dc.varID);
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
            
        for (FileSet fs : files) {
            for (DataColumn dc : fs.dataDefs) {
                lineOut = new StringBuilder();
                lineOut.append(dc.varName).append("\t");
                lineOut.append(dc.varID).append("\t");
                lineOut.append(dc.varDesc).append("\t");
                lineOut.append(dc.varUnit == null ? "." : dc.varUnit).append("\t");
                for (Entry<String, String> ent : dc.varValueDescriptions.entrySet()) {
                    lineOut.append(ent.getKey()).append("=").append(ent.getValue()).append(";");
                }
                if (dc.comment != null && !"".equals(dc.comment)) {
                    lineOut.append("\t").append(dc.comment);
                }
                writer.println(lineOut.toString());
            }
        }
        writer.flush();
        writer.close();
    }
    
    public static void main(String[] args) {
        String fileDir1 = "/home/pankrat2/shared/ARIC/phenos/dbGaP/50859/PhenoGenotypeFiles/RootStudyConsentSet_phs000280.ARIC_RootStudy.v3.p1.c1.HMB-IRB/PhenotypeFiles/";
        String fileDir2 = "/home/pankrat2/shared/ARIC/phenos/dbGaP/50865/PhenoGenotypeFiles/RootStudyConsentSet_phs000280.ARIC_RootStudy.v3.p1.c2.DS-CVD-IRB/PhenotypeFiles/";
        (new DBGAPMerge()).run(new String[]{fileDir1, fileDir2}, "/scratch.global/coleb/mergeOut.xln.gz", "/scratch.global/coleb/mergeMap.xln");
    }
    
}
