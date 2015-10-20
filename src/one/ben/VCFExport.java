package one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import common.Files;
import common.ext;


public class VCFExport {
    
    private VCFExport() {}
    
    static void exportToXLN(String vcfFile) throws IOException {
        
        // read file up til '#CHROM'
        //    collect all ##INFO lines, parse out into set, split pipes [prefix with orig. ##INFO id?]
        // write header
        // read lines
        
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
                    String[] desc = description.substring(description.indexOf("'") + 1).split("[\\s]*[\\(\\|][\\s]*"); 
                    multiTags.put(id, desc);
                    for (String descTag : desc) {
                        infoTags.add(descTag); 
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
            infoParts = vcfParts[7].split(";");
            for (String infoTag : infoParts) {
                String[] infoTagParts = infoTag.split("=");
                if (infoTagParts.length == 1) {
                    outLine.append("\t").append(infoTagParts[0]);
                } else if (multiTags.containsKey(infoTagParts[0])) {
                    String[] tagParts = infoTagParts[1].split("[\\(\\|]"); 
                    for (String tagPart : tagParts) {
                        outLine.append("\t").append("".equals(tagPart) ? "." : tagPart);
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

}
