package one.ben;

import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;

public class SNPEffTesting {
    
    public static void main(String[] args1) {
        String file = "D:\\All_20150605.vcf.gz";
        String config = "C:\\Users\\cole0482\\Downloads\\snpEff_latest_core\\snpEff\\snpEff.config";
        String[] args = {"ann", "-v", "-c", config, "hg19", file};
        SnpEff snpEff = new SnpEff(args);
        snpEff.run();
    }
    
    
}
