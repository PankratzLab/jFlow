package one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import common.Files;
import common.ext;
import cnv.filesys.Pedigree;
import cnv.filesys.Project;
import cnv.filesys.Pedigree.PedigreeEntry;
import cnv.qc.MarkerMetrics;
import cnv.qc.MendelErrors;
import cnv.qc.MendelErrors.MendelErrorCheck;
import cnv.qc.SampleQC;
import cnv.var.SampleData;

public class PLINKMendelianChecker {
    
    static class GenomeLoader {
        
        HashMap<String, HashMap<String, String>> pairData = new HashMap<String, HashMap<String,String>>();
        
        private GenomeLoader() {
            pairData = new HashMap<String, HashMap<String,String>>();
        }
        
        static GenomeLoader run(String genomeFile/*, HashMap<String, String[]> pedToFAMO*/, ArrayList<String[]> pairs, HashMap<String, String> dnaLookup) {
            BufferedReader reader;
            String line;
            String[] temp;
            
            GenomeLoader gl = new GenomeLoader();
            
//            boolean[] found = Array.booleanArray(pairs.size(), false);
//            int foundCount = 0;
            
            HashSet<String> dnas = new HashSet<String>();
            for (int i = 0; i < pairs.size(); i++) {
                String[] pair = pairs.get(i);
                String dna1 = dnaLookup.get(pair[0]);
                String dna2 = dnaLookup.get(pair[1]);
                dnas.add(dna1);
                dnas.add(dna2);
            }
            
            try {
                reader = Files.getAppropriateReader(genomeFile);
                line = reader.readLine().trim(); // read header
                temp = null;
                while ((line = reader.readLine()) != null) {
                    temp = line.split("[\\s]+");
//                temp[1] // FID1
//                temp[2] // IID1
//                temp[3] // FID2
//                temp[4] // IID2
//                temp[5] // RT    
//                temp[6] // EZ      
//                temp[7] // Z0      
//                temp[8] // Z1      
//                temp[9] // Z2  
//                temp[10] // PI_HAT 
//                temp[11] // PHE       
//                temp[12] // DST     
//                temp[13] // PPC   
//                temp[14] // RATIO
                    int FID1Index = 1;
                    int FID2Index = 3;
                    
                    if (dnas.contains(temp[FID1Index]) && dnas.contains(temp[FID2Index])) {
                        HashMap<String, String> data = gl.pairData.get(temp[FID1Index]);
                        if (data == null) {
                            data = new HashMap<String, String>();
                            gl.pairData.put(temp[FID1Index], data);
                        }
                        data.put(temp[FID2Index], line);
                        
                        HashMap<String, String> data2 = gl.pairData.get(temp[FID2Index]);
                        if (data2 == null) {
                            data2 = new HashMap<String, String>();
                            gl.pairData.put(temp[FID2Index], data2);
                        }
                        data2.put(temp[FID1Index], line);
                        
//                        for (int i = 0; i < pairs.size(); i++) {
//                            if (found[i]) continue;
//                            String[] pair = pairs.get(i);
//                            String dna1 = dnaLookup.get(pair[0]);
//                            String dna2 = dnaLookup.get(pair[1]);
//                            
//                            if ((temp[FID1Index].equals(dna1) && temp[FID2Index].equals(dna2)) || (temp[FID1Index].equals(dna2) && temp[FID2Index].equals(dna1))) {
//                                HashMap<String, String> data = gl.pairData.get(temp[FID1Index]);
//                                if (data == null) {
//                                    data = new HashMap<String, String>();
//                                    gl.pairData.put(temp[FID1Index], data);
//                                }
//                                data.put(temp[FID2Index], line);
//                                found[i] = true;
//                                foundCount++;
//                                if (pairs.size() < 100 || (pairs.size() > 100 && pairs.size() < 1000 && foundCount % 100 == 0) || (pairs.size() > 2000 && foundCount % 500 == 0)) {
//                                    System.out.println(ext.getTime() + "]\t Found " + foundCount + " of " + pairs.size() + " data sets");
//                                }
//                                break;
//                            }
//                        }
//                        if (Array.booleanArraySum(found) == found.length) {
//                            break;
//                        }
                    }
                }
                reader.close();
                reader = null;
            } catch (IOException e) {
                e.printStackTrace();
            }
            
            return gl;
        }
        
    }
    
    static class MendelLoader {
        HashMap<String, ArrayList<String>> errorMarkersMapFather;
        HashMap<String, ArrayList<String>> errorMarkersMapMother;
        private MendelLoader() {
            this.errorMarkersMapFather = new HashMap<String, ArrayList<String>>();
            this.errorMarkersMapMother = new HashMap<String, ArrayList<String>>();
        }
        
        static MendelLoader run(String mendelFile) {
            MendelLoader ml = new MendelLoader();
            BufferedReader reader;
            String line;
            String[] temp;
            try {
                reader = Files.getAppropriateReader(mendelFile);
                line = reader.readLine();
                while ((line = reader.readLine()) != null) {
                    temp = line.split("[\\s]+");
                    
                    MendelErrorCheck mec = new MendelErrors((byte)(Integer.parseInt(temp[1])), 
                                                                    -1, 
                                                                    (byte)(Integer.parseInt(temp[8])), 
                                                                    (byte)(Integer.parseInt(temp[9])), 
                                                                    (byte)(Integer.parseInt(temp[10]))).checkMendelError();
//                    MarkerName  
//                    Chr 
//                    Position    
//                    FID 
//                    IID 
//                    DNA
//                    FA_DNA  
//                    MO_DNA  
//                    AB  
//                    FA_AB   
//                    MO_AB
                    
                    if (mec.hasFaMendelError()) {
                        ArrayList<String> mkrs = ml.errorMarkersMapFather.get(temp[5]);
                        if (mkrs == null) {
                            mkrs = new ArrayList<String>();
                            ml.errorMarkersMapFather.put(temp[5], mkrs);
                        }
                        mkrs.add(temp[0]);
                    }
                    if (mec.hasMoMendelError()) {
                        ArrayList<String> mkrs = ml.errorMarkersMapMother.get(temp[5]);
                        if (mkrs == null) {
                            mkrs = new ArrayList<String>();
                            ml.errorMarkersMapMother.put(temp[5], mkrs);
                        }
                        mkrs.add(temp[0]);
                    }
                    
                }
                reader.close();
                reader = null;
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            return ml;
        }
        
    }
    
    final Project project;
    final String pedFile;
    final String genomeFile; // from property - cluster.genome.gz
    final String mendelFile; // from MarkerMetrics.fullQC - outputs a property file with an appended name
    final String outDir;
    
    public PLINKMendelianChecker(Project project) {
        this.project = project;
        this.pedFile = project.PEDIGREE_FILENAME.getValue(false, false);
        this.mendelFile = ext.rootOf(project.MARKER_METRICS_FILENAME.getValue(true, false), false) + MarkerMetrics.DEFAULT_MENDEL_FILE_SUFFIX;
        this.genomeFile = project.GENOME_CLUSTER_FILENAME.getValue();
        this.outDir = project.PROJECT_DIRECTORY.getValue(false, false);
    }
    
    public PLINKMendelianChecker(String pedFile, String mendelFile, String genomeFile, String outDir) {
        this.project = null;
        this.pedFile = pedFile;
        this.mendelFile = mendelFile;
        this.genomeFile = genomeFile;
        this.outDir = outDir;
    }
    
    
    void run() {
        Pedigree ped;
        SampleData sampleData = null;
        SampleQC sampQC = null;
        String[] samples = null;
        HashSet<String> pedDNA;
        HashMap<String, String> dnaLookup, idLookup;
        HashMap<String, String[]> pedToFAMO;
        HashMap<String, ArrayList<String>> childrenMap;
        HashMap<String, Integer> qcIndexMap = null;
        PrintWriter writer;
        GenomeLoader gl = null;
        MendelLoader ml = null;
        
        ped = Pedigree.loadPedigree(this.project, this.pedFile);
        if (ped == null) {
            // TODO Error
            return;
        }
        
        pedDNA = new HashSet<String>();
        pedToFAMO = new HashMap<String, String[]>();
        childrenMap = new HashMap<String, ArrayList<String>>();
        dnaLookup = new HashMap<String, String>();
        idLookup = new HashMap<String, String>();
        
        for (PedigreeEntry pe : ped.getPedigreeEntries()) {
            if (pe == null) {
                continue;
            }
            dnaLookup.put(pe.getFID() + "\t" + pe.getIID(), pe.getiDNA());
            idLookup.put(pe.getiDNA(), pe.getFID() + "\t" + pe.getIID());
            pedDNA.add(pe.getiDNA());
            pedToFAMO.put(pe.getiDNA(), new String[]{"0".equals(pe.getFA()) ? "." : pe.getFID() + "\t" + pe.getFA(), "0".equals(pe.getMO()) ? "." : pe.getFID() + "\t" + pe.getMO()});
            if (!"0".equals(pe.getFA())) {
                ArrayList<String> children = childrenMap.get(pe.getFID() + "\t" + pe.getFA());
                if (children == null) {
                    children = new ArrayList<String>();
                    childrenMap.put(pe.getFID() + "\t" + pe.getFA(), children);
                }
                children.add(pe.getFID() + "\t" + pe.getIID());
            }
            if (!"0".equals(pe.getMO())) {
                ArrayList<String> children = childrenMap.get(pe.getFID() + "\t" + pe.getMO());
                if (children == null) {
                    children = new ArrayList<String>();
                    childrenMap.put(pe.getFID() + "\t" + pe.getMO(), children);
                }
                children.add(pe.getFID() + "\t" + pe.getIID());
            }
        }
        
        ArrayList<String[]> pairs = new ArrayList<String[]>();
        for (java.util.Map.Entry<String, ArrayList<String>> childrenList : childrenMap.entrySet()) {
            String fidiid = childrenList.getKey();
            
            for (String childFIDIID : childrenList.getValue()) {
                pairs.add(new String[]{fidiid, childFIDIID});
                pairs.add(new String[]{childFIDIID, fidiid});

                String[] famo = pedToFAMO.get(dnaLookup.get(childFIDIID));
                ArrayList<String> spouseChildren = null;
                if (famo[0].equals(fidiid) && !".".equals(famo[1])) {
                    spouseChildren = childrenMap.get(famo[1]);
                } else if (famo[1].equals(fidiid) && !".".equals(famo[0])) {
                    spouseChildren = childrenMap.get(famo[0]);
                }
                
                HashSet<String> sameParentSibs = new HashSet<String>(childrenList.getValue());
                for (String otherChild : spouseChildren) {
                    if (otherChild.equals(childFIDIID)) continue;
                    if (sameParentSibs.contains(otherChild)) {
                        sameParentSibs.remove(otherChild);
                    }
                    pairs.add(new String[]{childFIDIID, otherChild});
                }
                for (String halfSib : sameParentSibs) {
                    if (halfSib.equals(childFIDIID)) continue;
                    pairs.add(new String[]{childFIDIID, halfSib});
                }
            }
        }
        
        if (genomeFile != null && (new File(genomeFile)).exists()) {
            System.out.println(ext.getTime() + "]\tLoading GenomicCluster data...");
            gl = GenomeLoader.run(genomeFile/*, pedToFAMO*/, pairs, dnaLookup);
        }
        if (mendelFile != null && (new File(mendelFile)).exists()) {
            System.out.println(ext.getTime() + "]\tLoading MendelianError data...");
            ml = MendelLoader.run(mendelFile);
        }
        
        if (project != null) {
            System.out.println(ext.getTime() + "]\tLoading Project data...");
            samples = project.getSamples();
            sampleData = project.getSampleData(0, false);
            sampQC = SampleQC.loadSampleQC(project);
            qcIndexMap = new HashMap<String, Integer>();
            if (sampQC != null) {
                for (int i = 0; i < sampQC.getSamples().length; i++) {
                    qcIndexMap.put(sampQC.getSamples()[i], i);
                }
            }
        }
        
        StringBuilder sb;
        sb = new StringBuilder()
            .append("FID").append("\t")
            .append("IID").append("\t")
            .append("FA_IID").append("\t")
            .append("MO_IID").append("\t")
            .append("SEX").append("\t")
            .append("DNA").append("\t")
            .append("FA_DNA").append("\t")
            .append("MO_DNA").append("\t");
        if (gl != null) {
            sb.append("IBD0_FATHER").append("\t")
                .append("IBD1_FATHER").append("\t")
                .append("IBD2_FATHER").append("\t")
                .append("PIHAT_FATHER").append("\t")
                .append("IBD0_MOTHER").append("\t")
                .append("IBD1_MOTHER").append("\t")
                .append("IBD2_MOTHER").append("\t")
                .append("PIHAT_MOTHER").append("\t");
        }
        if (ml != null) {
            sb.append("MendelianErrorsFather").append("\t")
                .append("MendelianErrorsMother").append("\t");
        }
        
        if (sampQC != null) {
            sb.append("LRR_SD_CHILD").append("\t")
                .append("LRR_SD_FATHER").append("\t")
                .append("LRR_SD_MOTHER").append("\t")
                .append("CALLRATE_CHILD").append("\t")
                .append("CALLRATE_FATHER").append("\t")
                .append("CALLRATE_MOTHER").append("\t");
        }
        
        if (sampleData != null) {
            sb.append("EXCLUDE_CHILD").append("\t")
                .append("EXCLUDE_FATHER").append("\t")
                .append("EXCLUDE_MOTHER").append("\t");
        }
        
        System.out.println(ext.getTime() + "]\tWriting result data...");
        
        String outputHeader = sb.toString();
        writer = Files.getAppropriateWriter(outDir + "trios.xln");
        writer.println(outputHeader);
        for (int i = 0; i < ped.getPedigreeEntries().length; i++) {
            PedigreeEntry pe = ped.getPedigreeEntries()[i];
            if (pe == null) continue;
            sb = new StringBuilder();
            sb.append(pe.getFID()).append("\t")
                .append(pe.getIID()).append("\t")
                .append(pe.getFA()).append("\t")
                .append(pe.getMO()).append("\t")
                .append(pe.getSEX()).append("\t")
                .append(pe.getiDNA()).append("\t");
            
            String faDNA;
            if (pe.getFaDNAIndex() == -1 || samples == null) {
                if ("0".equals(pe.getFA()) || dnaLookup.get(pe.getFID() + "\t" + pe.getFA()) == null) {
                    faDNA = ".";
                } else {
                    faDNA = dnaLookup.get(pe.getFID() + "\t" + pe.getFA());
                }
            } else {
                faDNA = samples[pe.getFaDNAIndex()];
            }
            sb.append(faDNA).append("\t");
            String moDNA;
            if (pe.getMoDNAIndex() == -1 || samples == null) {
                if ("0".equals(pe.getMO()) || dnaLookup.get(pe.getFID() + "\t" + pe.getMO()) == null) {
                    moDNA = ".";
                } else {
                    moDNA = dnaLookup.get(pe.getFID() + "\t" + pe.getMO());
                }
            } else {
                moDNA = samples[pe.getMoDNAIndex()];
            }
            sb.append(moDNA).append("\t");
            
            if (gl != null) {
                
                HashMap<String, String> genoLines = gl.pairData.get(pe.getiDNA());
                if (genoLines == null) {
                    System.out.println("No genomic data for " + pe.getiDNA());
                }
                if (genoLines == null || ".".equals(faDNA) || genoLines.get(faDNA) == null) {
                    sb.append(".").append("\t")
                        .append(".").append("\t")
                        .append(".").append("\t")
                        .append(".").append("\t");   
                } else {
                    String[] tmpGL = genoLines.get(faDNA).split("[\\s]+");
                    sb.append(tmpGL[7]).append("\t")
                        .append(tmpGL[8]).append("\t")
                        .append(tmpGL[9]).append("\t")
                        .append(tmpGL[10]).append("\t");
                }
                if (genoLines == null || ".".equals(moDNA) || genoLines.get(moDNA) == null) {
                    sb.append(".").append("\t")
                        .append(".").append("\t")
                        .append(".").append("\t")
                        .append(".").append("\t");   
                } else {
                    String[] tmpGL = genoLines.get(moDNA).split("[\\s]+");
                    sb.append(tmpGL[7]).append("\t")
                        .append(tmpGL[8]).append("\t")
                        .append(tmpGL[9]).append("\t")
                        .append(tmpGL[10]).append("\t");
                }
            }
            
            if (ml != null) {
                sb.append(ml.errorMarkersMapFather.get(pe.getiDNA()).size()).append("\t")
                    .append(ml.errorMarkersMapMother.get(pe.getiDNA()).size()).append("\t");
            }
            
            if (sampQC != null) {
                double[] data = sampQC.getDataFor("LRR_SD");
                Integer indexInt = qcIndexMap.get(pe.getiDNA());
                Integer faIndex = null;
                if (pe.getFaDNAIndex() == -1) {
                    if (dnaLookup.get(pe.getFID() + "\t" + pe.getFA()) != null) {
                        faIndex = qcIndexMap.get(dnaLookup.get(pe.getFID() + "\t" + pe.getFA()));
                    }
                } else {
                    faIndex = qcIndexMap.get(pe.getFaDNAIndex());
                }
                Integer moIndex = null;
                if (pe.getMoDNAIndex() == -1) {
                    if (dnaLookup.get(pe.getFID() + "\t" + pe.getMO()) != null) {
                        moIndex = qcIndexMap.get(dnaLookup.get(pe.getFID() + "\t" + pe.getMO()));
                    }
                } else {
                    moIndex = qcIndexMap.get(pe.getMoDNAIndex());
                }
                if (indexInt == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[indexInt.intValue()]).append("\t");
                }
                if (faIndex == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[faIndex.intValue()]).append("\t");
                }
                if (moIndex == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[moIndex.intValue()]).append("\t");
                }
                data = sampQC.getDataFor("Genotype_callrate");
                if (indexInt == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[indexInt.intValue()]).append("\t");
                }
                if (faIndex == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[faIndex.intValue()]).append("\t");
                }
                if (moIndex == null) {
                    sb.append(".").append("\t");
                } else {
                    sb.append(data[moIndex.intValue()]).append("\t");
                }
            }
            
            if (sampleData != null) {
                sb.append(sampleData.individualShouldBeExcluded(pe.getiDNA()) ? "1" : "0").append("\t");
                if (pe.getFaDNAIndex() == -1 || samples == null) {
                    if (dnaLookup.get(pe.getFID() + "\t" + pe.getFA()) != null) {
                        sb.append(sampleData.individualShouldBeExcluded(dnaLookup.get(pe.getFID() + "\t" + pe.getFA())) ? "1" : "0").append("\t");
                    } else {
                        sb.append(".").append("\t");
                    }
                } else {
                    sb.append(sampleData.individualShouldBeExcluded(samples[pe.getFaDNAIndex()]) ? "1" : "0").append("\t");
                }
                if (pe.getMoDNAIndex() == -1 || samples == null) {
                    if (dnaLookup.get(pe.getFID() + "\t" + pe.getMO()) != null) {
                        sb.append(sampleData.individualShouldBeExcluded(dnaLookup.get(pe.getFID() + "\t" + pe.getMO())) ? "1" : "0").append("\t");
                    } else {
                        sb.append(".").append("\t");
                    }
                } else {
                    sb.append(sampleData.individualShouldBeExcluded(samples[pe.getMoDNAIndex()]) ? "1" : "0").append("\t");
                }
            }
            
            writer.println(sb.toString());
        }
        
        writer.flush();
        writer.close();
        
        writeFamily(gl, pedToFAMO, childrenMap, dnaLookup, outDir);
    }
    
    private void writeFamily(GenomeLoader gl, HashMap<String, String[]> pedToFAMO, HashMap<String, ArrayList<String>> childrenMap, HashMap<String, String> dnaLookup, String outDir2) {
        PrintWriter writer;
        StringBuilder sb;
        
        writer = Files.getAppropriateWriter(outDir + "family.xln");
        
        String header = "FID1\tIID1\tFID2\tIID2\tPUTATIVE_REL\t";
        if (gl != null) {
            header += "IBD0" + "\t"
                      + "IBD1" + "\t"
                      + "IBD2" + "\t"
                      + "PIHAT";
        }
        writer.println(header);
        
        for (java.util.Map.Entry<String, ArrayList<String>> childrenList : childrenMap.entrySet()) {
            String fidiid = childrenList.getKey();
            
            for (String childFIDIID : childrenList.getValue()) {
                
                sb = new StringBuilder();
                sb.append(fidiid).append("\t");
                sb.append(childFIDIID).append("\t")
                    .append("PO").append("\t");
                if (gl != null) {
                    String dna = dnaLookup.get(fidiid);
                    if (dna == null) {
                        sb.append(".\t.\t.\t.\t");
                    } else {
                        HashMap<String, String> genoData = gl.pairData.get(dna);
                        String dna2 = dnaLookup.get(childFIDIID);
                        if (genoData != null && dna2 != null) {
                            String genoDataLine = genoData.get(dna2);
                            if (genoDataLine != null) {
                                String[] tmpGL = genoDataLine.split("[\\s]+");
                                sb.append(tmpGL[7]).append("\t")
                                    .append(tmpGL[8]).append("\t")
                                    .append(tmpGL[9]).append("\t")
                                    .append(tmpGL[10]).append("\t");
                            } else {
                                sb.append(".\t.\t.\t.\t");
                            }
                        } else {
                            sb.append(".\t.\t.\t.\t");
                        }
                    }
                }
                writer.println(sb.toString());
                
                String[] famo = pedToFAMO.get(dnaLookup.get(childFIDIID));
                ArrayList<String> spouseChildren = null;
                if (famo[0].equals(fidiid) && !".".equals(famo[1])) {
                    spouseChildren = childrenMap.get(famo[1]);
                } else if (famo[1].equals(fidiid) && !".".equals(famo[0])) {
                    spouseChildren = childrenMap.get(famo[0]);
                }
                
                HashSet<String> sameParentSibs = new HashSet<String>(childrenList.getValue());
                for (String otherChild : spouseChildren) {
                    if (otherChild.equals(childFIDIID)) continue;
                    sb = new StringBuilder();
                    sb.append(childFIDIID).append("\t")
                        .append(otherChild).append("\t");
                    if (sameParentSibs.contains(otherChild)) {
                        sameParentSibs.remove(otherChild);
                        sb.append("SIB").append("\t");
                    } else {
                        sb.append("HALFSIB").append("\t");
                    }
                    if (gl != null) {
                        String dna = dnaLookup.get(fidiid);
                        if (dna == null) {
                            sb.append(".\t.\t.\t.\t");
                        } else {
                            HashMap<String, String> genoData = gl.pairData.get(dna);
                            String dna2 = dnaLookup.get(childFIDIID);
                            if (genoData != null && dna2 != null) {
                                String genoDataLine = genoData.get(dna2);
                                if (genoDataLine != null) {
                                    String[] tmpGL = genoDataLine.split("[\\s]+");
                                    sb.append(tmpGL[7]).append("\t")
                                        .append(tmpGL[8]).append("\t")
                                        .append(tmpGL[9]).append("\t")
                                        .append(tmpGL[10]).append("\t");
                                } else {
                                    sb.append(".\t.\t.\t.\t");
                                }
                            } else {
                                sb.append(".\t.\t.\t.\t");
                            }
                        }
                    }
                }
                for (String halfSib : sameParentSibs) {
                    if (halfSib.equals(childFIDIID)) continue;
                    sb = new StringBuilder();
                    sb.append(childFIDIID).append("\t")
                        .append(halfSib).append("\t")
                        .append("HALFSIB").append("\t");
                    if (gl != null) {
                        String dna = dnaLookup.get(fidiid);
                        if (dna == null) {
                            sb.append(".\t.\t.\t.\t");
                        } else {
                            HashMap<String, String> genoData = gl.pairData.get(dna);
                            String dna2 = dnaLookup.get(childFIDIID);
                            if (genoData != null && dna2 != null) {
                                String genoDataLine = genoData.get(dna2);
                                if (genoDataLine != null) {
                                    String[] tmpGL = genoDataLine.split("[\\s]+");
                                    sb.append(tmpGL[7]).append("\t")
                                        .append(tmpGL[8]).append("\t")
                                        .append(tmpGL[9]).append("\t")
                                        .append(tmpGL[10]).append("\t");
                                } else {
                                    sb.append(".\t.\t.\t.\t");
                                }
                            } else {
                                sb.append(".\t.\t.\t.\t");
                            }
                        }
                    }
                }
                
            }
            
        }
        writer.flush();
        writer.close();
        
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String projFile = null;
        String ped = null;
        String mendel = null;
        String genomic = null;
        String out = null;

        String usage = "\n" + 
                       "one.ben.PLINKMendelianChecker requires 0-1 arguments\n" + 
                       "   (1) Project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
                       "  OR \n" + 
                       "   (1) File with pedigree data (i.e. pedigree=pedigree.dat (not the default))\n" + 
                       "   (2) File with Mendelian Error data (i.e. mendel=markerQualityChecks.mendel (not the default))\n" + 
                       "   (2) File with genomic cluster data (i.e. genomic=cluster.genome.gz (not the default))\n" + 
                       "   (3) Directory of output (i.e. out=/path/to/dir/ (not the default))\n"+
                       "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("proj=")) {
                projFile = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("pedigree=")) {
                ped = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("mendel=")) {
                mendel = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("genomic=")) {
                genomic = args[i].split("=")[1];
                numArgs--;
            } else if (args[i].startsWith("out=")) {
                out = args[i].split("=")[1];
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
            if (projFile != null) {
                (new PLINKMendelianChecker(new Project(projFile, false))).run();
            } else if (ped != null && mendel != null && genomic != null && out != null) {
                (new PLINKMendelianChecker(ped, mendel, genomic, out)).run();
            } else {
                System.err.println(usage);
                System.exit(1);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
