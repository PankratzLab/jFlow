package one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import common.Files;
import common.ext;


public class VCFExport {
    
    private static final String TAG_SPLIT_SET = "[\\(\\|]";

    private VCFExport() {}
    
    public static void exportToXLN(String vcfFile) throws IOException {
        
        BufferedReader reader;
        PrintWriter writer;
        ArrayList<String> infoTags;
        HashMap<String, String[]> multiTags;
        StringBuilder header;
        String[] infoParts, vcfParts;
        
        reader = Files.getAppropriateReader(vcfFile);
        writer = Files.getAppropriateWriter(ext.rootOf(vcfFile, false) + ".xln");
        infoTags = new ArrayList<String>();
        multiTags = new HashMap<String, String[]>();
        String line = null;
        while ((line = reader.readLine()) != null && !line.startsWith("#CHROM")) {
            if (line.startsWith("##INFO")) {
                int ind1 = line.indexOf("ID=") + "ID=".length();
                int ind2 = line.indexOf(",", ind1);
                String id = line.substring(ind1, ind2);
                infoTags.add(id);
                ind1 = line.indexOf("Description=") + "Description=".length() + 1;
                ind2 = line.indexOf("\"", ind1);
                String description = line.substring(ind1, ind2);
                if (description.contains("|")) {
                    infoTags.remove(id);
                    String[] desc = description.substring(description.indexOf("'") + 1).split("[\\s]*" + TAG_SPLIT_SET + "[\\s]*"); 
                    multiTags.put(id, desc);
                    for (String descTag : desc) {
                        infoTags.add(id + "_" + descTag); 
                    }
                }
            }
        }
        if (line == null) {
            // error 
            return;
        }
        header = new StringBuilder("rsID\tChr\tPos\tRef\tAlt\tQual\tFilter");
        for (String tag : infoTags) {
            header.append("\t").append(tag.replaceAll("[\\[\\(\\]\\)\\']", "").trim());
        }
        writer.println(header);
        
        while((line = reader.readLine()) != null) {
            vcfParts = line.split("\t");
            StringBuilder outLine = new StringBuilder();
            outLine.append(vcfParts[2]).append("\t") // rs
                    .append(vcfParts[0]).append("\t") // chr
                    .append(vcfParts[1]).append("\t") // pos
                    .append(vcfParts[3]).append("\t") // ref
                    .append(vcfParts[4]).append("\t") // alt
                    .append(vcfParts[5]).append("\t") // qual
                    .append(vcfParts[6]); // filt
            if (vcfParts[7].startsWith(".;")) {
                vcfParts[7] = vcfParts[7].substring(2);
            }
            infoParts = vcfParts[7].split(";");
            for (String infoTag : infoParts) {
                String[] infoTagParts = infoTag.split("=");
                if (infoTagParts.length == 1) {
                    outLine.append("\t").append(infoTagParts[0]);
                } else if (multiTags.containsKey(infoTagParts[0])) {
                    String[] tagParts = infoTagParts[1].split(TAG_SPLIT_SET, -1); 
                    for (String tagPart : tagParts) {
                        tagPart = tagPart.trim();
                        if ("".equals(tagPart)) {
                            tagPart = ".";
                        } else if (tagPart.contains("(") && !tagPart.contains(")")) {
                            tagPart = tagPart.replaceAll("\\(", "");
                        } else if (tagPart.contains(")") && !tagPart.contains("(")) {
                            tagPart = tagPart.replaceAll("\\)", "");
                        } else if (tagPart.contains("[") && !tagPart.contains("[")) {
                            tagPart = tagPart.replaceAll("\\[", "");
                        } else if (tagPart.contains("]") && !tagPart.contains("]")) {
                            tagPart = tagPart.replaceAll("\\]", "");
                        }
                        
                        outLine.append("\t").append(tagPart);
                    }
                } else {
                    outLine.append("\t").append(infoTagParts[1]);
                }
            }
            writer.println(outLine);
            outLine = null;
            
        }
        writer.flush();
        writer.close();
        reader.close();
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = "F:/variants.vcf";//"VCFExport.dat";

        String usage = "\n" + 
                        "one.ben.VCFExport requires 0-1 arguments\n" + 
                        "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("file=")) {
                filename = args[i].split("=")[1];
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
            exportToXLN(filename);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
