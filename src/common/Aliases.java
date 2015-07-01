package common;

public class Aliases {
	public static final String[] MARKER_NAMES = {"MarkerName", "Marker", "Name", "SNP", "SNPID", "SNP.id", "rs_id", "Variant", "AnalysisUnit", "Variant_ID", "SingleVariant", "BurdenTests", "RSID" };
	public static final String[] GENE_UNITS = {"Gene", "SKATgene"};
	
	public static final String[] CHRS = {"Chr", "Chromosome", "CHROM"};
	public static final String[] POSITIONS = {"Position", "position", "pos", "BP", "MapInfo", "PositionOfFirstMarkerInGene"};
	public static final String[] CENTIMORGANS = {"centiMorgans", "cM"};

	public static final String[][] ALLELES = {
		{"coded_all", "A1", "Al1", "Allele1", "ALT", "Effect_allele", "EA"},
		{"noncoded_all", "A2", "Al2", "Allele2", "REF", "OTHER", "Reference_allele", "NEA", "Other_allele"},		
	};
	public static final String[] EFFECTS = {"beta", "beta_SNP_add", "Effect"};
	public static final String[] STD_ERRS = {"se", "StdErr", "sebeta_SNP_add"};
	
	public static final String[] PVALUES = {"pval", "P", "p-val", "p-value", "Pvalue", "mbpval", "minPval"};
	public static final String[] NS = {"N", "NMISS", "sampleN", "n_total"};
	
	public static final String[] ALLELE_FREQS = {"freq", "AlleleFreq", "A1Freq", "AF", "AAF", "MAF", "sampleMAF", "Effect_allele_frequency", "EAF"};
	
	public static final String[] IMPUTATION_EFFICIENCY = {"imp_info", "rsq"};
}
