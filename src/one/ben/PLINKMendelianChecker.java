package one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import common.Files;
import common.ext;
import cnv.filesys.Pedigree;
import cnv.filesys.Project;
import cnv.filesys.Pedigree.PedigreeEntry;
import cnv.qc.SampleQC;
import cnv.var.SampleData;


public class PLINKMendelianChecker {
    
    String genomeFile; // from property - cluster.genome.gz
    String mendelFile; // from MarkerMetrics.fullQC - outputs a property file with an appended name
    // TODO double check why we need/want the mendel file - only extra info is specific markers with genotype disparities; perhaps a listing of these for trios?
    Project project; // may be null? - no sample data, all files must be specified by hand, NO DNA LIST
    
    public PLINKMendelianChecker(Project project) {
        this.project = project;
        mendelFile = ext.rootOf(this.project.MARKER_METRICS_FILENAME.getValue(true, false), false) + "_mendel.error";
        
    }
    
    void run() {
        Pedigree ped;
        SampleData sampleData;
        SampleQC sampQC;
        String line;
        String[] samples, temp;
        HashSet<String> pedDNA;
        HashMap<String, String[]> pedToFAMO, genomeLines;
        BufferedReader reader;
        
        
        ped = project.loadPedigree();
        if (ped == null) {
            // TODO Error
            return;
        }
       
        if (project != null) {
            sampleData = project.getSampleData(0, false);
            sampQC = SampleQC.loadSampleQC(project);
        }
        
        samples = project.getSamples();
        pedDNA = new HashSet<String>();
        pedToFAMO = new HashMap<String, String[]>();
        genomeLines = new HashMap<String, String[]>();
        
        for (PedigreeEntry pe : ped.getPedigreeEntries()) {
            if (!pe.getiDNA().equals(samples[pe.getiDNAIndex()])) {
                // TODO ERROR/DNA MISMATCH
            }
            pedDNA.add(pe.getiDNA());
            pedToFAMO.put(pe.getiDNA(), new String[]{samples[pe.getFaDNAIndex()], samples[pe.getMoDNAIndex()]});
            genomeLines.put(pe.getiDNA(), new String[2]);
        }
        
        try {
            reader = Files.getAppropriateReader(genomeFile);
            line = reader.readLine(); // read header
            temp = null;
            while ((line = reader.readLine()) != null) {
                temp = line.split("[\\s]+");
                
                if (pedDNA.contains(temp[0])) {
                    if (pedToFAMO.get(temp[0])[0].equals(temp[2])) {
                        // Father
                        genomeLines.get(temp[0])[0] = line;
                    } else if (pedToFAMO.get(temp[0])[1].equals(temp[2])) {
                        // Mother
                        genomeLines.get(temp[0])[1] = line;
                    } else {
                        continue; // neither, discard
                    }
                }
//            temp[0] // FID1
//            temp[1] // IID1
//            temp[2] // FID2
//            temp[3] // IID2
//            temp[4] // RT    
//            temp[5] // EZ      
//            temp[6] // Z0      
//            temp[7] // Z1      
//            temp[8] // Z2  
//            temp[9] // PI_HAT 
//            temp[10] // PHE       
//            temp[11] // DST     
//            temp[12] // PPC   
//            temp[13] // RATIO
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        // read mendel file, same thing
        
        
        
    }
    
    
}
