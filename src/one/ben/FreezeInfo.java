package one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import common.Array;
import common.Files;
import common.ext;

public class FreezeInfo {
    
//    "CHROM","POS","REF","ALT","func_region","SNP","SKATgene","sc_nonsynSplice","sc_functional","MAF","AAF"
//    "10",100154981,"T","C","ncRNA_exonic","10:100154981","MIR1287",FALSE,FALSE,NA,NA
//    "10",100154995,"G","A","ncRNA_exonic","10:100154995","MIR1287",FALSE,FALSE,NA,NA
    
    
//"SNP","CHROM","POS","REF","ALT","ANNOVAR_ucsc_precedent_consequence","ANNOVAR_ucsc_precedent_gene","INDEL","sc_exonic","sc_lof","sc_indel","sc_nonsynSplice","sc_damaging","sc_indel_coding","sc_functional","inFreeze5","SKATgene","dmg_total","sc_damaging_orig"
//"10:100003242:T:G","10",100003242,"T","G","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003242:T:G",0,FALSE
//"10:100003302:G:A","10",100003302,"G","A","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003302:G:A",0,FALSE
//"10:100003304:A:G","10",100003304,"A","G","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003304:A:G",0,FALSE
//"10:100003747:C:G","10",100003747,"C","G","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003747:C:G",0,FALSE
//"10:100003753:T:C","10",100003753,"T","C","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003753:T:C",0,FALSE
//"10:100003766:G:C","10",100003766,"G","C","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003766:G:C",0,FALSE
//"10:100003785:T:C","10",100003785,"T","C","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003785:T:C",0,FALSE
//"10:100003796:T:A","10",100003796,"T","A","intronic","R3HCC1L",NA,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,"10:100003796:T:A",0,FALSE

    static String FRZ_4_SRC = "F:/freeze4/snpinfo_ChargeSFreeze3Freeze4_ESP_RS_ERF_Broad_Analytic_04112014.csv";
    static String FRZ_5_SRC = "F:/freeze4/snpinfo_WES_v13_Analytic_10202015.csv.gz";
    
    static String ALTERED_OUT = "F:/freeze4/freeze4.altered.out";
    static String ALTERED_ALL_OUT = "F:/freeze4/freeze4.altered.all.out";
    static String NO_MATCH_OUT = "F:/freeze4/freeze4.noMatch.out";
    static String NO_MATCH_ALL_OUT = "F:/freeze4/freeze4.noMatch.all.out";
    
    static String LOG_OUT = "F:/freeze4/freeze4.log";
    static String LOG_ALL_OUT = "F:/freeze4/freeze4.all.log";
    
    static String HEADER = "CHROM\tPOS\tREF\tALT\tfunc_region\tSNP\tSKATgene\tsc_nonsynSplice\tsc_functional\tMAF\tAAF";
    
    static void run(boolean ALL) throws IOException {
        HashMap<String, ArrayList<String[]>> frz5 = load();
    
        BufferedReader frz4reader = Files.getAppropriateReader(FRZ_4_SRC);
        PrintWriter finalWriter = Files.getAppropriateWriter(ALL ? ALTERED_ALL_OUT : ALTERED_OUT);
        PrintWriter noMatchWriter = Files.getAppropriateWriter(ALL ? NO_MATCH_ALL_OUT : NO_MATCH_OUT);
    
        String line = null;
        String[] parts;
    
        frz4reader.readLine();
        finalWriter.println(HEADER);
        
        int countExonicToExonic = 0;
        int countNonToExonic = 0;
        int countFuncToFunc = 0;
        int countFuncToNonFunc = 0;
        int countNonFuncToFunc = 0;
        int countNonFuncToNonFunc = 0;
        int missing = 0;
        int noMatch = 0;
        
        while ((line = frz4reader.readLine()) != null) {
            parts = ext.splitCommasIntelligently(line, true, null);
            if (!ALL && !line.contains("exonic")) {
                finalWriter.println(Array.toStr(parts, "\t"));
                continue;
            }
            int ind = parts[5].indexOf(":", 4);
            String mkr = ind == -1 ? parts[5] : parts[5].substring(0, ind);
            String ref4 = parts[2];
            String alt4 = parts[3];
            String cons4 = parts[4];
            String gene4 = parts[6];
            String func4 = parts[8];
            ArrayList<String[]> frz5Info = frz5.get(mkr);
            
            String[] frz5data = null;
            boolean exact = false;
            if (frz5Info != null) {
                if (!frz5Info.isEmpty()) {
                    for (String[] data : frz5Info) {
                        String gene5 = data[1];
                        if (gene4.equalsIgnoreCase(gene5)) {
                            frz5data = data;
                            exact = true;
                            break;
                        }
                    }
                    if (frz5data == null) {
                        for (String[] data : frz5Info) {
                            String ref5 = data[4];
                            String alt5 = data[5];
                            String gene5 = data[1];
                            if (gene5.contains(gene4) || gene4.contains(gene5)) {
                                if (ref4.equals(ref5) && alt4.equals(alt5)) {
                                    frz5data = data;
                                    break;
                                }
                            } else if (gene4.split("-")[0].equals(gene5.split("-")[0])) {
                                if (ref4.equals(ref5) && alt4.equals(alt5)) {
                                    frz5data = data;
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    System.out.println("ERROR - empty"); // TODO
                }
            } else {
                System.out.println("ERROR - missing"); // TODO
                missing++;
            }
            if (frz5data != null) {
                String cons5 = frz5data[0];
                String func5 = frz5data[3];
                
                boolean ex4 = cons4.contains("exonic");
                boolean ex5 = cons5.contains("exonic");
                boolean f4 = func4.equalsIgnoreCase("true");
                boolean f5 = func5.equalsIgnoreCase("true");
                
                if (ex4 && ex5) {
                    countExonicToExonic++;
                } else if (!ex4 && ex5) {
                    countNonToExonic++;
                } else if (f4 && f5) {
                    countFuncToFunc++;
                } else if (f4 && !f5) {
                    countFuncToNonFunc++;
                } else if (!f4 && f5) {
                    countNonFuncToFunc++;
                } else if (!f4 &&!f5) {
                    countNonFuncToNonFunc++;
                }
                
                if (exact) {
                    parts[8] = func5;
                    parts[4] = cons5;
                }
                
            } else {
                noMatch++;
                noMatchWriter.println(line);
                if (frz5Info != null) {
                    for (String[] data : frz5Info) {
                        noMatchWriter.println(Array.toStr(data, "\t"));
                    }
                }
                noMatchWriter.println();
                noMatchWriter.println();
            }
            
            finalWriter.println(Array.toStr(parts, "\t"));
        }
        finalWriter.flush();
        finalWriter.close();
        noMatchWriter.flush();
        noMatchWriter.close();
        frz4reader.close();        
        
        PrintWriter writer = Files.getAppropriateWriter(ALL ? LOG_ALL_OUT : LOG_OUT);
        writer.println("Exonic -> Exonic : " + countExonicToExonic);
        writer.println("NonExonic -> Exonic : " + countNonToExonic);
        writer.println("Functional -> Functional : " + countFuncToFunc);
        writer.println("Functional -> NonFunctional : " + countFuncToNonFunc);
        writer.println("NonFunctional -> Functional : " + countNonFuncToFunc);
        writer.println("NonFunctional -> NonFunctional : " + countNonFuncToNonFunc);
        writer.println("Missing --> " + missing);
        writer.println("No Match --> " + noMatch);
        writer.flush();
        writer.close();
    }
    
    static HashMap<String, ArrayList<String[]>> load() throws IOException {
        HashMap<String, ArrayList<String[]>> frz5Data = new HashMap<String, ArrayList<String[]>>();
        
        BufferedReader frz5reader = Files.getAppropriateReader(FRZ_5_SRC);
        String line = null;
        String[] parts;
        frz5reader.readLine();
        while ((line = frz5reader.readLine()) != null) {
            parts = ext.splitCommasIntelligently(line, true, null);
            int ind = parts[0].indexOf(":", 4);
            String mkr = ind == -1 ? parts[0] : parts[0].substring(0, ind);
            ArrayList<String[]> data = frz5Data.get(mkr);
            if (data == null) {
                data = new ArrayList<String[]>();
                frz5Data.put(mkr, data);
            }
            data.add(new String[]{parts[5], parts[6], parts[8], parts[14], parts[3], parts[4]});
        }
        
        return frz5Data;
    }
    
    public static void main(String[] args) {
        try {
            run(false);
            System.out.println("2");
            run(true);
            System.out.println("complete");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    
}
