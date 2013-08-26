package filesys;

import gwas.Metal;

public class MetaAnalysisParams {
	public static final String[] STUDIES = {"ARIC", "CHS", "FHS", "ESP6800", "RS", "Rergo3"};
	public static final String[] STUDY_GROUPS = {"CHARGE", "CHARGE", "CHARGE", "ESP", "RS"}; // +ESP
	
	public static final String[][][] ALL_PHENOTYPES = {
		{
			{"Fibrinogen", "fibrinogen", ".lnFB.", ".ln_fib."}, 
			{"F7", ".FVII.", "_FVII_"}, 
			{"F8", ".FVIII.", "_FVIII_"}, 
			{"vWF", "VWF"}
		},
		{
			{"CRP"}, 
		},
		
	};

	public static final String[][] RACES = {
		{"Whites", ".EA.", "_EA_", ".EU.", "_EU_"}, 
		{"Blacks", ".AA.", "_AA_"}, 
	};

	public static final int[][] DEFAULT_SAMPLE_SIZES = {
		{985, 628, 255, 2011}, // Fibr
		{965, 630, 248, 1221}, // F7
		{984, 626, 0, 1296}, // F8
		{985, 0, 249, 1184}, // vWF
	};
	public static final int[][] FREEZE3_SAMPLE_SIZES = {
		{3168, 741, 499, 2013}, // Fibr
		{3092, 744, 416, 1222}, // F7
		{3166, 736, 0, 1296}, // F8
		{3168, 0, 416, 1185}, // vWF
	};

	public static final String SNP_INFO_FILE = "snpinfo_ChargeSFreeze3_ESP_05212013_slimmest.RData";
	public static final String SNP_NAMES = "SNP";
	public static final String CHROM_NAME = "CHROM";
	public static final String GENE_NAME = "SKATgene";

	public static final String[][] METHODS = {{"SingleSNP", "singleSnpRes", ".LR."}, {"T5Count", ".T5."}, {"T5MB", "MBT5"}};
	public static final String[][] UNIT_OF_ANALYSIS = {Metal.MARKER_NAMES, Metal.GENE_UNITS, Metal.GENE_UNITS};
	public static final boolean[] SINGLE_VARIANTS = {true, false, false};
	public static final boolean[] WEIGHTED = {true, false, false};
	public static final String[] GROUPS = {"SingleVariant", "BurdenTests", "BurdenTests"};	

	public static final String[][] GROUP_ANNOTATION_PARAMS = {
		{},
//		{"../SNPInfo_ExomeFreeze2_120810_aafSlim.csv 'SNP' 'gene' 'AAF'=CHARGE_AF"},
		{},
	};
	
}
