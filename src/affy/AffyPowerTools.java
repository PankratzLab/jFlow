package affy;
import java.io.*;
import java.util.Vector;
import common.Array;
import common.Files;
import common.CmdLine;
import common.ext;

// target sketch? 
public class AffyPowerTools {
    
    /**apt-probeset-genotype -c lib/GenomeWideSNP_6.cdf -a birdseed --read-models-birdseed lib/GenomeWideSNP_6.birdseed.models 
    * --special-snps lib/GenomeWideSNP_6.specialSNPs --out-dir apt --cel-files listfile
    *call rate >96%
    *SD of normalized intensity (LRR) <0.35
    *the GC wave factor of LRR was between −0.02 and 0.02 were accepted
    *CNV call count <80 
    *
    */
    
    public static void main(String[] args) {
        System.out.println("affy");
        String affyLib = "/project/scratch/normDat/affy/lib/";
        String outDir = "/project/scratch/normDat/output/";
        String pennCNVbin = "/project/scratch/normDat/pennCNV/bin/";
        String pennCNVlib = "/project/scratch/normDat/pennCNV/lib/";
        //String intitialCelList = "/project/scratch/normDat/lists/bigListFinal.txt";
        String finalCelList = "/project/scratch/normDat/lists/bigListFinal.txt";
        String pbsOutDir = "/project/scratch/normDat/";
        
        String affyChunk = "cel";
        String affyChunkProbe = "probes" ;
        String affyChunkProbeAll = "allProbes" ;
        int numJobs = 2;
        String celLists = "/project/scratch/normDat/lists/";
        
        String affyDetectPBS= pbsOutDir + "Detect.pbs";
        String affyQCPBS= pbsOutDir + "QC.pbs";
        String affyKColPBS= pbsOutDir + "KCol.pbs";
        String affyGenoPBS= pbsOutDir + "GT.pbs";
        String affySumPBS= pbsOutDir + "Sum.pbs";
        String affyDetctPBS= pbsOutDir + "Detect.pbs";
        String affyPennGenoClustPBS= pbsOutDir + "PenGClust.pbs";
        String affyPennGenoLRRBAFPBS= pbsOutDir + "PennLRR_BAF.pbs";
        String affyCDF = affyLib + "GenomeWideSNP_6.cdf";
        String affyBirdseedModel = affyLib + "GenomeWideSNP_6.birdseed-v2.models";
        String affySpecialSnps = affyLib + "GenomeWideSNP_6.specialSNPs";
        String affyHapMapQuant = affyLib + "hapmap.quant-norm.normalization-target.txt";
        String affyRefInput = affyLib + "GenomeWideSNP_6.hapmap270.na33.r1.a5.ref";
        String affyChrX = affyLib + "GenomeWideSNP_6.chrXprobes";
        String affyChrY = affyLib + "GenomeWideSNP_6.chrYprobes";
        String affyQCA = affyLib + "GenomeWideSNP_6.r2.qca";
        String affyQCC = affyLib + "GenomeWideSNP_6.r2.qcc";
        String annoFile = affyLib + "GenomeWideSNP_6.na33.annot.db";
        String affyGenoQC = affyLib + "apt-geno-qc";
        String affyGenotype = affyLib + "apt-probeset-genotype";
        String affySummarize = affyLib + "apt-probeset-summarize";
        String qcOut =  outDir + "qcOut.txt";
        String sexFile = outDir + "file_sex.txt"; 
        String birdseedReport =  "birdseed-v2.report.txt"; 
        String birdseedCalls =  "birdseed-v2.calls.txt"; 
        String birdseedConf =  "birdseed-v2.confidences.txt"; 
        String quantNorm =  "sort-all-quant-norm.pm-only.med-polish.expr.summary.txt";
//        
//        String affyPennCnv = pennCNV + "/gw6/bin/";
//        String affyPennCnvLib = pennCNV + "gw6/lib/";
        String locFile = pennCNVlib + "affygw6.hg18.pfb";
        String genoCluster =  "gw6.genocluster"; 
        String genoLrrBaf =  "gw6.lrr_baf.txt"; 
        String pennSplitOut = "signalListFile.txt";
        String pennHmm = pennCNVlib + "affygw6.hmm";
        String detect_cnv = pennCNVbin + "penncnv/";
        
        
        //QC params
        double callRateCut = .95;
        
        //works
//        String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " + 
//                affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY  + " --cel-files " 
//                + intitialCelList + " --out-file " + qcOut;    
//        String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " + 
//        affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY  + " --out-file $wdir/qcOut.qc" + " *.CEL";    
//        //works
//        String genotypeCommand = affyGenotype + " --summaries --write-models -c " + affyCDF + " -a birdseed-v2 --read-models-birdseed " +
//                affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir " +
//                outDir + " --cel-files " + intitialCelList + " --chrX-probes " + affyChrX +" --chrY-probes " + 
//                affyChrY ;
        //works
//        String summarizeCommand =  "apt-probeset-summarize --cdf-file " + affyCDF +
//                " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "+
//                affyHapMapQuant + " --out-dir " + outDir + " --cel-files " + finalCelList;
        //works
//        String copyNumberWorkflow = "apt-copynumber-workflow -v 1 --cdf-file " + affyCDF +
//                " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY +  " --special-snps " 
//                + affySpecialSnps + " --annotation-file " + annoFile + " --reference-input " +
//                affyRefInput + " -out-dir " + outDir + " --cel-files " + finalCelList + " --text-output true"; 
        //works
//        String generate_affy_geno_cluster  = "perl " + affyPennCnv+"generate_affy_geno_cluster.pl "+
//                birdseedCalls + " " + birdseedConf + " " + quantNorm + " -locfile " +
//                locFile + " -sexfile  " + sexFile + " -out " + genoCluster; 
        //works
//        String lrrBafCalc = "perl " + affyPennCnv+"normalize_affy_geno_cluster.pl "+
//                genoCluster + " " + quantNorm + " -locfile " + locFile +
//                " -out " + genoLrrBaf;
//        //works
//        String splitSignalFile = "perl " + pennCNV+"kcolumn.pl "+
//                genoLrrBaf + " split 2 -tab -head 3 --name_by_header -out gw6 --filenameout " +pennSplitOut ;
//        //works        
//        String detectCNV = "perl " + pennCNV+"detect_cnv.pl --test -hmm " + pennHmm +
//                " -pfb " + locFile + " -list " + pennSplitOut + " -log gw6.log -out gw6.rawcnv";
//                
                
        //System.out.println(qcCommand);
//        CmdLine.run(qcCommand , outDir );
//        
        //genFinalCelList(qcOut , finalCelList , callRateCut);
        //System.out.println(genotypeCommand);
//        CmdLine.run(genotypeCommand , outDir);
//        
        //System.out.println(summarizeCommand);
//        CmdLine.run(summarizeCommand , outDir);
//        
        //genSexFile(birdseedReport, sexFile);
//        //System.out.println(copyNumberWorkflow);
////        CmdLine.run(copyNumberWorkflow , outDir);
//        System.out.println(generate_affy_geno_cluster);
//        CmdLine.run(generate_affy_geno_cluster , outDir);
        
        //System.out.println(lrrBafCalc);
        //CmdLine.run(lrrBafCalc , outDir);
        //run command from output directory
        
        //System.out.println(splitSignalFile);
        //CmdLine.run(splitSignalFile , outDir);
        
        //System.out.println(detectCNV);
        //CmdLine.run(detectCNV , outDir);
        
        String[] affyQCJobs = new String[numJobs];
        String[] affyGenoJobs = new String[numJobs];
        String[] affySumJobs = new String[numJobs];
        String[] affyGenoClusterJobs = new String[numJobs];
        String[] LRRBAFJobs = new String[numJobs];
        String[] kColumn = new String[numJobs];
        String[] detectCNVs = new String[numJobs];
        
        
        for (int i = 0; i<numJobs; i++) {
            String batchNum = "" + i;
            if(i<10){
                batchNum = "0" + i;
            }
            
            
            String qc = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " + 
                affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY  + " --out-file "+outDir + "qcOut_" + batchNum 
                + ".qc --cel-files " + celLists + affyChunk +batchNum;
            affyQCJobs[i] = qc;
            
            String genotypeCommand = affyGenotype + " -c " + affyCDF + " -a birdseed-v2 --set-gender-method cn-probe-chrXY-ratio --read-models-birdseed " +
                    affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir " +
                    outDir + "genoTypeOut" + batchNum + "/ --cel-files " + finalCelList + " --chrX-probes " + affyChrX +" --chrY-probes " + 
                    affyChrY + " --probeset-ids " + celLists + affyChunkProbe +batchNum ;
            affyGenoJobs[i] = genotypeCommand;
            
            String summarizeCommand =  affySummarize +" --cdf-file " + affyCDF +
                    " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "+
                    affyHapMapQuant + " --out-dir " + outDir + "summarOut" + batchNum +
                    "/ --cel-files " + finalCelList
                    + " --probeset-ids " + celLists + affyChunkProbeAll +batchNum;
            
            
            affySumJobs[i] = summarizeCommand;
            
            // Include sex file later
            String generate_affy_geno_cluster  = "perl " + pennCNVbin+"generate_affy_geno_cluster.pl "+
                    outDir + "genoTypeOut" + batchNum +"/" +  birdseedCalls 
                    + " " + outDir + "genoTypeOut" + batchNum +"/" + birdseedConf + " " 
                    + outDir +quantNorm + " -locfile " +
                    pennCNVlib + locFile +  " -out " + outDir + "genoTypeOut" + batchNum + "/" + genoCluster; 
            
            affyGenoClusterJobs[i] = generate_affy_geno_cluster;
            
            String lrrBafCalc = "perl " + pennCNVbin+"normalize_affy_geno_cluster.pl "+
                    outDir + "genoTypeOut" + batchNum + "/" + genoCluster + " " + outDir +quantNorm 
                    + " -locfile " + pennCNVlib + locFile + 
                    " -out " + outDir + "genoTypeOut" + batchNum + "/" +genoLrrBaf;
            
            LRRBAFJobs[i] = lrrBafCalc;
            
            String splitSignalFile = "perl " + pennCNVbin + "kcolumn.pl "+
            		outDir + "genoTypeOut" + batchNum + "/" + genoLrrBaf + 
            		" split 2 -tab -head 3 --name -out " + outDir + "Kcol/Kcol" + 
            		batchNum +"/gw6_split --filenameout " +
            		outDir + "Kcol/Kcol" + batchNum + "/" +pennSplitOut ;
            kColumn[i] = splitSignalFile;
            
            String detectCNV = "perl " + detect_cnv +"detect_cnv.pl --test -hmm " + pennHmm +
                  " -pfb " + locFile + " -log " +outDir + "detectCNV" +
                  batchNum + "/gw6.log -out " +outDir + "detectCNV" + batchNum + "/gw6.rawcnv" +
                  " -list " + celLists + "indDetect" +batchNum;;
                  
          
            detectCNVs[i] = detectCNV;
            //System.out.println(qc);
            //System.out.println(generate_affy_geno_cluster);
            System.out.println(detectCNV);
        }
        
        //qsubMultiple(affyQCPBS , affyQCJobs, 3000 , 4);
        Files.qsubMultiple(affyGenoPBS , affyGenoJobs, 3000 , 4);
        Files.qsubMultiple(affySumPBS , affySumJobs, 3000 , 4);
        Files.qsubMultiple(affyPennGenoClustPBS , affyGenoClusterJobs, 3000 , 4);
        Files.qsubMultiple(affyPennGenoLRRBAFPBS , LRRBAFJobs, 3000 , 4);
        Files.qsubMultiple(affyKColPBS , kColumn, 3000 , 4);
        Files.qsubMultiple(affyDetectPBS , detectCNVs, 3000 , 4);
        
        
        
    }
    
    
    
   
    private static void genFinalCelList (String aQCReport , String finalCelList , double callRate){
        PrintWriter writer;
        BufferedReader reader;
        String[] line;
        try {
            writer = new PrintWriter(new FileWriter(finalCelList));
            reader = new BufferedReader(new FileReader(aQCReport));
            while (reader.ready()) {
                line = reader.readLine().split("\t");
                    if (line[0].matches("#.*|cel_files")){
                            continue;    
                    }
                    else if(Double.parseDouble(line[1])  > callRate ) {
                        writer.write(line[0] + "\n");
                    }
                    else{
                        System.out.println(line[1]);
                    }                    
            }
            reader.close();
            writer.close();
        } 
        catch (IOException ioe) {
            System.err.println("Error reading file \"" + aQCReport + "\"");
            System.exit(2);
        }
    }
        
        
    
    
    

    private static String genSexFile(String birdseedConfFile, String aSexFile ) {
        PrintWriter writer;
        BufferedReader reader;
        String[] line;
        try {
            writer = new PrintWriter(new FileWriter(aSexFile));
            reader = new BufferedReader(new FileReader(birdseedConfFile));
            while (reader.ready()) {
                line = reader.readLine().split("\t");
                    if (line[0].matches("#.*|cel_files")){
                            continue;    
                    }
                    else if (line[1].matches("unknown")){
                        continue;
                    }
                    else if(line[1].equals("female")) {
                        writer.write(line[0] + "\t0\n");
                        
                    }
                    else if(line[1].equals("male")) {
                        writer.write(line[0] + "\t1\n");
                    }
                    else{
                        System.out.println(line[1]);
                    }                    
            }
            reader.close();
            writer.close();
        } 
        catch (IOException ioe) {
            System.err.println("Error reading file \"" + birdseedConfFile + "\"");
            System.exit(2);
        }
        return aSexFile;
    }
}



