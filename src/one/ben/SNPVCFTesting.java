package one.ben;

import filesys.SerialHash;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import seq.manage.ReferenceGenome;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import cnv.annotation.AnnotationParser;
import cnv.annotation.MarkerAnnotationLoader;
import cnv.annotation.AnnotationFileLoader.QUERY_ORDER;
import cnv.filesys.Project;

public class SNPVCFTesting {

    static String vcfFile = "D:/All_20150605.vcf.gz";
    static String unzippedVCFFile = "D:/All_20150605.vcf";
    static String trimmedVCFFile = "D:/All_20150605_trimmed.vcf.gz";
    static String numberedVCFFile = "D:/All_20150605_noRS.vcf.gz";
    static String trimmedVCFFile2 = "D:/All_20150605_trimmed_rerun.vcf.gz";
    static String numberedVCFFile2 = "D:/All_20150605_noRS_rerun.vcf.gz";

    static String mergedFile = "D:/RsMerge.ser";
    static String mergedOut = "D:/RsMerge.vcf.gz";
    static String mergedOut2 = "D:/RsMerge_rerun.vcf.gz";
    
    static String UNIX_unzippedVCFFile = "/scratch.global/coleb/DBSNP/All_20150605.vcf";
    static String UNIX_numberedVCFFile = "/scratch.global/coleb/DBSNP/All_20150605_noRS.vcf.gz";
    static String UNIX_numberedVCFFile2 = "/scratch.global/coleb/DBSNP/All_20150605_noRS_rerun.vcf.gz";
//    static String tbiFile = "D:/All_20150605.vcf.gz.tbi";
    
    static String papu = "D:/All_20150605_papu.vcf";
    static String papu_out = "D:/All_20150605_papu_trimmedNoRS.vcf.gz";
    static String papu_out2 = "D:/All_20150605_papu_trimmedNoRS_rerun.vcf.gz";
    
    static String[] header = new String[]{
        "##fileformat=VCFv4.0",
        "##fileDate=20150605",
        "##source=dbSNP",
        "##dbSNP_BUILD_ID=144",
        "##reference=GRCh37.p13",
        "##phasing=partial",
        "##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)\">",
        "##INFO=<ID=PM,Number=0,Type=Flag,Description=\"Variant is Precious(Clinical,Pubmed Cited)\">",
        "##INFO=<ID=NSF,Number=0,Type=Flag,Description=\"Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44\">",
        "##INFO=<ID=NSM,Number=0,Type=Flag,Description=\"Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42\">",
        "##INFO=<ID=NSN,Number=0,Type=Flag,Description=\"Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41\">",
        "##INFO=<ID=REF,Number=0,Type=Flag,Description=\"Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8\">",
        "##INFO=<ID=SYN,Number=0,Type=Flag,Description=\"Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3\">",
        "##INFO=<ID=U3,Number=0,Type=Flag,Description=\"In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53\">",
        "##INFO=<ID=U5,Number=0,Type=Flag,Description=\"In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55\">",
        "##INFO=<ID=ASS,Number=0,Type=Flag,Description=\"In acceptor splice site FxnCode = 73\">",
        "##INFO=<ID=DSS,Number=0,Type=Flag,Description=\"In donor splice-site FxnCode = 75\">",
        "##INFO=<ID=INT,Number=0,Type=Flag,Description=\"In Intron FxnCode = 6\">",
        "##INFO=<ID=R3,Number=0,Type=Flag,Description=\"In 3' gene region FxnCode = 13\">",
        "##INFO=<ID=R5,Number=0,Type=Flag,Description=\"In 5' gene region FxnCode = 15\">",
        "##INFO=<ID=MUT,Number=0,Type=Flag,Description=\"Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources\">",
        "##INFO=<ID=CAF,Number=.,Type=String,Description=\"An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    };
    
    static String chrPosInfo = "##INFO=<ID=CHRPOS,Number=1,Type=String,Description=\"Chromosome and Position, separated by a colon (:)\">";
    
    static HashSet<String> validTags = new HashSet<String>(){{
            add("GENEINFO");
            add("PM");
            add("NSF");
            add("NSM");
            add("NSN");
            add("REF");
            add("SYN");
            add("U3");
            add("U5");
            add("ASS");
            add("DSS");
            add("INT");
            add("R3");
            add("R5");
            add("MUT");
            add("CAF");
        }
    };
            
    static void createTBI() {
        String refGenome = "D:/hg19_canonical.fa";
        SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(refGenome, new Logger()).getIndexedFastaSequenceFile().getSequenceDictionary();
        
//        VCFFileReader fileReader1 = new VCFFileReader(new File(trimmedVCFFile), false);
//        VariantContextWriter vcw = VCFOps.initWriter(trimmedVCFFile2, VCFOps.DEFUALT_WRITER_OPTIONS, samSequenceDictionary);
//        VCFHeader hdr = fileReader1.getFileHeader();
//        hdr.setSequenceDictionary(samSequenceDictionary);
//        vcw.writeHeader(hdr);
//        for (VariantContext vc : fileReader1) {
//            VariantContextBuilder vcB = new VariantContextBuilder(vc);
//            
//            String chr = vc.getChr();
//            if (ext.isValidInteger(chr)) {
//                chr = Positions.getChromosomeUCSC(Integer.parseInt(chr), true);
//            } else {
//                if (chr.equals("MT")) {
//                    chr = "chrM";
//                } else {
//                    chr = "chr" + chr;
//                }
//            }
//            
//            vcB.chr(chr);
//            vcw.add(vcB.make());
//        }
//        vcw.close();
//        fileReader1.close();
//
//        System.out.println("one");
        
        samSequenceDictionary = (new SAMSequenceDictionary());
        samSequenceDictionary.addSequence(new SAMSequenceRecord("chr1", (512 * 1024 * 1024)));
//        samSequenceDictionary.addSequence(new SAMSequenceRecord("chr2", (512 * 1024 * 1024)));

        VCFFileReader fileReader2 = new VCFFileReader(new File(mergedOut), false);
        VariantContextWriter vcw2 = VCFOps.initWriter(mergedOut2, VCFOps.DEFUALT_WRITER_OPTIONS, samSequenceDictionary);
        VCFHeader hdr = fileReader2.getFileHeader();
        
        hdr.setSequenceDictionary(samSequenceDictionary);
        vcw2.writeHeader(hdr);
        for (VariantContext vc : fileReader2) {
            
            VariantContextBuilder vcB = new VariantContextBuilder(vc);
            String chr = vc.getContig();
            if (ext.isValidInteger(chr)) {
                chr = Positions.getChromosomeUCSC(Integer.parseInt(chr), true);
            } else {
                if (chr.equals("MT")) {
                    chr = "chrM";
                } else {
                    chr = "chr" + chr;
                }
            }
            vcB.chr(chr);
            VariantContext vc2 = vcB.make();
            try {
                vcw2.add(vc2);
            } catch (Exception e) {
                e.printStackTrace();
                fileReader2.close();
                return;
            }
        }
        vcw2.close();
        fileReader2.close();
        
        System.out.println("done");
    }
    
    
    static void trimFileOnce() throws IOException {
        BufferedReader reader = Files.getAppropriateReader(vcfFile);
        PrintWriter writerTrimmed = Files.getAppropriateWriter(trimmedVCFFile);
        
        for (int i = 0; i < header.length; i++) {
            writerTrimmed.println(header[i]);
        }
        
        String line = null;
        boolean afterHeader = false;
        while((line = reader.readLine()) != null) {
            line = line.trim();
            if (!afterHeader) {
                if (line.startsWith("#CHROM")) {
                    afterHeader = true;
                }
                continue;
            }
            String[] lineParts = line.split("\t");
            String[] infoParts = lineParts[lineParts.length - 1].split(";");

            StringBuilder newLineTrimmed = new StringBuilder();
            for (int i = 0; i < lineParts.length - 1; i++) {
                newLineTrimmed.append(lineParts[i]).append("\t");
            }
            int tags = 0;
            for (String str : infoParts) {
                if (validTags.contains(str) || (str.split("=").length > 1 && validTags.contains(str.split("=")[0]))) {
                    newLineTrimmed.append(str).append(";");
                    tags++;
                }
            }
            if (newLineTrimmed.charAt(newLineTrimmed.length() - 1) == ';') {
                newLineTrimmed.deleteCharAt(newLineTrimmed.length() - 1);
            }
            if (tags == 0) {
                newLineTrimmed.append(".");
            }
            writerTrimmed.println(newLineTrimmed.toString());

        }
        writerTrimmed.flush();
        writerTrimmed.close();
        reader.close();
    }
    
    static void writeByRSNumber() throws IOException {
        System.out.println(ext.getTime() + "]\tCounting lines...");
        int lines = Files.countLines(papu, 0);
        System.out.println(ext.getTime() + "]\t" + lines + " lines in VCF file.");
        System.out.println(ext.getTime() + "]\tReading header...");
        BufferedReader reader = Files.getAppropriateReader(papu);
        boolean afterHeader = false;
        String line = null;
        int headerLineCount = 0;
        while((line = reader.readLine()) != null && !afterHeader) {
            line = line.trim();
            if (!afterHeader) {
                headerLineCount++;
                if (line.startsWith("#CHROM")) {
                    afterHeader = true;
                    break;
                }
                continue;
            }
        }
        System.out.println(ext.getTime() + "]\t" + headerLineCount + " header lines in VCF file.");
        
        int MAX_CONTIG_LEN = 512 * 1024 * 1024;
        
        int[] rsNumbers = new int[lines - headerLineCount];
        HashMap<Integer, String> rsNumberLineMap = new HashMap<Integer, String>();
        
        System.out.println(ext.getTime() + "]\tReading data...");
        int cnt = 0;
        HashSet<Integer> contigs = new HashSet<Integer>();
        while((line = reader.readLine()) != null) {
            line = line.trim();
            String[] lineParts = line.split("\t");
            String[] infoParts = lineParts[lineParts.length - 1].split(";");
            String rsStr = null;
            for (String k : infoParts) {
                if (k.startsWith("RS=")) {
                    rsStr = k.split("=")[1];
                }
            }
            if (rsStr == null) {
                System.err.println("Error - not found");
                continue;
            }
            int rsNumber = Integer.parseInt(rsStr);
            rsNumbers[cnt] = rsNumber;
            cnt++;
            StringBuilder newLineNonRS = new StringBuilder();
            String chrPos = "CHRPOS=" + lineParts[0] + ":" + lineParts[1];
            int contig = (rsNumber / MAX_CONTIG_LEN) + 1;
            if (!contigs.contains(contig)) {
                contigs.add(contig);
            }
            newLineNonRS.append(contig).append("\t");
            newLineNonRS.append(rsNumber % MAX_CONTIG_LEN).append("\t");
            newLineNonRS.append("rs" + rsNumber).append("\t");
            for (int i = 3; i < lineParts.length - 1; i++) {
                newLineNonRS.append(lineParts[i]).append("\t");
            }

            newLineNonRS.append(chrPos);
            for (String str : infoParts) {
                if (validTags.contains(str) || (str.split("=").length > 1 && validTags.contains(str.split("=")[0]))) {
                    newLineNonRS.append(";").append(str);
                }
            }
            rsNumberLineMap.put(rsNumber, newLineNonRS.toString());
        }
        System.out.println("Count: " + cnt);
        System.out.println("Length: " + rsNumbers.length);
        System.out.println("Contigs: " + contigs);
        
        System.out.println(ext.getTime() + "]\tSorting RS numbers...");
        Arrays.sort(rsNumbers);
        
        System.out.println(ext.getTime() + "]\tWriting new file...");
        PrintWriter writerNonRS = Files.getAppropriateWriter(papu_out);
        
        for (int i = 0; i < header.length; i++) {
            if (i == header.length - 1) {
                writerNonRS.println(chrPosInfo);
            }
            writerNonRS.println(header[i]);
        }
        for (int number : rsNumbers) {
            Integer key = Integer.valueOf(number);
            writerNonRS.println(rsNumberLineMap.get(key));
        }
        writerNonRS.flush();
        writerNonRS.close();
        
        System.out.println(ext.getTime() + "]\tComplete!"); 
    }
   
    static String hdr0 = "##fileformat=VCFv4.2";
    static String hdr1 = "##INFO=<ID=RSMRG,Number=1,Type=String,Description=\"RS number merged with this one\">";
    static String hdr2 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    
    static void processMerged() {
        Hashtable<Integer,Integer> mergeHash = SerialHash.loadSerializedIntHash(mergedFile);

        int MAX_CONTIG_LEN = 512 * 1024 * 1024;
        HashMap<Integer, String> rsNumberLineMap = new HashMap<Integer, String>();
        HashSet<Integer> contigs = new HashSet<Integer>();
        int[] rsNumbers = new int[mergeHash.size()];
        int index = 0;
        for (java.util.Map.Entry<Integer, Integer> mergeEntry : mergeHash.entrySet()) {
            Integer intKey = mergeEntry.getKey();
            rsNumbers[index++] = intKey;
            StringBuilder newLineNonRS = new StringBuilder();
            int contig = (intKey / MAX_CONTIG_LEN) + 1;
            if (!contigs.contains(contig)) {
                contigs.add(contig);
            }
            newLineNonRS.append(contig).append("\t");
            newLineNonRS.append(intKey % MAX_CONTIG_LEN).append("\t");
            newLineNonRS.append("rs" + intKey).append("\t");
            newLineNonRS.append("A\tG\t");
            newLineNonRS.append(".\t.\t");
            newLineNonRS.append("RSMRG=").append(mergeEntry.getValue());
            rsNumberLineMap.put(intKey, newLineNonRS.toString());
        }
        System.out.println("Length: " + rsNumbers.length);
        System.out.println("Contigs: " + contigs);
        
        Arrays.sort(rsNumbers);
        
        System.out.println(ext.getTime() + "]\tWriting new file...");
        PrintWriter writerNonRS = Files.getAppropriateWriter(mergedOut);
        
        writerNonRS.println(hdr0);
        writerNonRS.println(hdr1);
        writerNonRS.println(hdr2);
        for (int number : rsNumbers) {
            Integer key = Integer.valueOf(number);
            writerNonRS.println(rsNumberLineMap.get(key));
        }
        writerNonRS.flush();
        writerNonRS.close();
        
        System.out.println(ext.getTime() + "]\tComplete!"); 
    }
    
    
    static String[] testMarkers = new String[]{
        "rs12175489",
        "rs9469003",
        "rs7080373",
        "rs2524279",
        "rs2247296",
        "rs1503859",
        "rs5911809",
        "rs1470383",
        "rs13427710",
        "rs6735220",
        "rs2965174",
        "rs609418",
        "rs4465523",
        "rs13379232",
        "rs6994361",
        "rs1413710",
        "rs4934434",
        "rs933150",
        "rs10158752",
        "rs9803750",
        "rs4660123",
        "rs4372063",
        "rs1641661",
        "rs754118",
        "rs6020624",
    };
    
    static void run() {
        VCFFileReader reader = new VCFFileReader(numberedVCFFile2, true);
        
        for (String marker : testMarkers) {
            // skip check for RS
            int number = Integer.parseInt(marker.substring(2));
            
            int chrom = number / (512 * 1024 * 1024) + 1;
            System.out.println(marker);
            CloseableIterator<VariantContext> vcIter = reader.query("chr" + chrom, number, number);
            while(vcIter.hasNext()) {
                VariantContext vc = vcIter.next();
                Object val = vc.getAttribute("CHRPOS");
                System.out.println("\t" + vc.getID());
            }
            
        }
        
    }

    
    public static void main(String[] args) {
        processMerged();
//        run();
//        try {
//            writeByRSNumber();
            createTBI();
//            trimFileOnce();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }
    
}
