package one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import common.Array;
import common.Files;
import common.ext;

public class CARDIAConcatenator {

//    String fileDir1 = ext.verifyDirFormat("F:/CARDIA processing/consent 1/");
//    String fileDir2 = ext.verifyDirFormat("F:/CARDIA processing/consent 2/");
//
//    String outDir = ext.verifyDirFormat("F:/CARDIA processing/all/");
    String fileDir1 = ext.verifyDirFormat("/home/pankarne/pankratz/CARDIA/CARDIA.whites.1000G.imputed/imputed/consent_level_1/");
    String fileDir2 = ext.verifyDirFormat("/home/pankarne/pankratz/CARDIA/CARDIA.whites.1000G.imputed/imputed/consent_level_2/");
  
    String outDir = ext.verifyDirFormat("/home/pankarne/pankratz/CARDIA/CARDIA.whites.1000G.imputed/imputed/consent_level_all/");

    String fileIDs = "ids_chr#.txt";

    String fileFormats[] = new String[] {"out_chr#comb.cardia_plus_filtered_chr#.dose.gz" };

    public void run() {
        for (String formatTemplate : fileFormats) {
            for (int chrInd = 1; chrInd < 23; chrInd++) {
                System.out.println(ext.getTime() + "]\tProcessing Chr" + chrInd);
                String format = formatTemplate.replaceAll("#", chrInd + "");
                try {
                    BufferedReader reader1 = Files.getAppropriateReader(fileDir1 + format);
                    BufferedReader reader2 = Files.getAppropriateReader(fileDir2 + format);
                    PrintWriter writer = Files.getAppropriateWriter(outDir + format);
                    PrintWriter writerIDs = Files.getAppropriateWriter(outDir + fileIDs.replace("#", chrInd + ""));

                    String line1 = reader1.readLine();
                    String line2 = reader2.readLine();
                    
                    writer.println(combineHeader(line1, line2, writerIDs));

                    while ((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null) {
                        writer.println(combineLine(line1, line2));
                        // TODO non-matching lines .... in-memory hash instead of line-by-line?
                        // TODO test for nonmatching first - combineLine will throw an exception if mismatched
                    }

                    if (line1 != null) {
                        PrintWriter excess1 = Files.getAppropriateWriter(outDir + "chr" + chrInd + "_excess_c1.out");
                        excess1.println(line1);
                        while((line1 = reader1.readLine()) != null) {
                            excess1.println(line1);
                        }
                        excess1.flush();
                        excess1.close();
//                        throw new RuntimeException("Extra content in file 1");
                    }

                    if (line2 != null) {
                        PrintWriter excess2 = Files.getAppropriateWriter(outDir + "chr" + chrInd + "_excess_c2.out");
                        excess2.println(line2);
                        while((line2 = reader2.readLine()) != null) {
                            excess2.println(line2);
                        }
                        excess2.flush();
                        excess2.close();
//                        throw new RuntimeException("Extra content in file 2");
                    }
                    
                    writer.flush();
                    writerIDs.flush();
                    
                    writer.close();
                    writerIDs.close();
                    
                    writer = null;
                    writerIDs = null;
                    
                    reader1.close();
                    reader2.close();
                    
                    System.out.println(ext.getTime() + "]\tFinished Processing Chr" + chrInd);
                } catch (FileNotFoundException e) {
                    // e.printStackTrace();
                    System.out.println(e.getMessage());
                } catch (IOException e) {
                    // e.printStackTrace();
                    System.out.println(e.getMessage());
                }
            }
        }

    }

    private String combineLine(String line1, String line2) throws RuntimeException {
        String[] line1Parts = line1.split("[\\s]+");
        String[] line2Parts = line2.split("[\\s]+");
        
        if (!line1Parts[0].trim().equals(line2Parts[0].trim())) {
            throw new RuntimeException("Error - mismatched marker names: <"+line1Parts[0]+"> and <"+line2Parts[0]+">");
        }
        if (!line1Parts[1].trim().equalsIgnoreCase(line2Parts[1].trim())) {
            throw new RuntimeException("Error - mismatched allele 1: <"+line1Parts[1]+"> and <"+line2Parts[1]+">");
        }
        if (!line1Parts[2].trim().equalsIgnoreCase(line2Parts[2].trim())) {
            throw new RuntimeException("Error - mismatched allele 2: <"+line1Parts[2]+"> and <"+line2Parts[2]+">");
        }
        
        
        String[] newLine = new String[line1Parts.length + line2Parts.length - 3];
        newLine[0] = line1Parts[0];
        newLine[1] = line1Parts[1];
        newLine[2] = line1Parts[2];
        
        for (int i = 3; i < line1Parts.length; i++) {
            newLine[i] = line1Parts[i];
        }
        for (int i = 3; i < line2Parts.length; i++) {
            newLine[line1Parts.length + (i - 3)] = line2Parts[i];
        }
        
        return Array.toStr(newLine, " ");
    }

    private String combineHeader(String line1, String line2, PrintWriter writerIDs) {
        String[] line1Parts = line1.split("[\\s]+");
        String[] line2Parts = line2.split("[\\s]+");

        String[] newLine = new String[line1Parts.length + line2Parts.length - 3];
        newLine[0] = line1Parts[0];
        newLine[1] = line1Parts[1];
        newLine[2] = line1Parts[2];

        for (int i = 3; i < line1Parts.length; i++) {
            newLine[i] = line1Parts[i];
            writerIDs.println(line1Parts[i]);
        }
        for (int i = 3; i < line2Parts.length; i++) {
            newLine[line1Parts.length + (i - 3)] = line2Parts[i];
            writerIDs.println(line2Parts[i]);
        }
        
        return Array.toStr(newLine, " ");
    }
    
    public static void main(String[] args) {
        (new CARDIAConcatenator()).run();
    }
    
}
