package common;

public class Positions {
	public static final String[] CHR_CODES = {"Un", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "XY", "MT"};
	public static final int[] CENTROMERE_MIDPOINTS = {-1, 124300000, 93300000, 91700000, 50700000, 47700000, 60500000, 59100000, 45200000, 51800000, 40300000, 52900000, 35400000, 16000000, 15600000, 17000000, 38200000, 22700000, 16100000, 28500000, 27100000, 12300000, 11800000, 60000000, 11500000, -1, -1};
	public static final int[] CHROMOSOME_LENGTHS = {-1, 247200000, 242800000, 199500000, 191300000, 180900000, 170900000, 158900000, 146300000, 140300000, 135400000, 134500000, 132300000, 114200000, 106400000, 100400000, 88900000, 78700000, 76200000, 63900000, 62500000, 47000000, 49600000, 154600000, 57800000, 154900000, 20000, -1};
//	public static final String[] CENTROMERE_MIDPOINT_SEGS = {"chr1:125,000,000", "chr2:93,000,000", "chr3:91,000,000", "chr4:50,000,000", "chr5:48,000,000", "chr6:60,000,000", "chr7:59,500,000", "chr8:45,000,000", "chr9:53,000,000", "chr10:40,000,000", "chr11:54,000,000", "chr12:35,000,000", "chr13:15,600,000", "chr14:16,000,000", "chr15:16,000,000", "chr16:37,500,000", "chr17:22,700,000", "chr18:16,500,000", "chr19:28,300,000", "chr20:27,000,000", "chr21:11,500,000", "chr22:13,000,000", "chrX:61,000,000", "chrY:11,800,000"};
	public static final String[] CENTROMERE_MIDPOINT_SEGS = {"chr1:124,300,000", "chr2:93,300,000", "chr3:91,700,000", "chr4:50,700,000", "chr5:47,700,000", "chr6:60,500,000", "chr7:59,100,000", "chr8:45,200,000", "chr9:51,800,000", "chr10:40,300,000", "chr11:52,900,000", "chr12:35,400,000", "chr13:16,000,000", "chr14:15,600,000", "chr15:17,000,000", "chr16:38,200,000", "chr17:22,700,000", "chr18:16,100,000", "chr19:28,500,000", "chr20:27,100,000", "chr21:12,300,000", "chr22:11,800,000", "chrX:59,500,000", "chrY:11,800,000"};
	public static final int[][] CENTROMERE_BOUNDARIES_FROM_SNPS = {{-1, -1}, // determined using the Omni map
																{121186705, 141477078},
															    {91682828, 94691119},
															    {90576572, 94994003},
															    {49353840, 52355341},
															    {46440100, 49442477},
															    {58887738, 61940508},
															    {58023003, 61060840},
															    {43956838, 46958273},
															    {46992793, 65231255},
															    {39194226, 41677980},
															    {51447829, 54468566},
															    {34747647, 36144166},
															    {10384996, 17922259},
															    {0, 18070422},
															    {0, 18275409},
															    {35141900, 44943958},
															    {22699477, 22707804},
															    {15399965, 16770020},
															    {24421125, 32424385},
															    {26267039, 28039018},
															    {10208131, 13260157},
															    {0, 14430353},
															    {58588152, 61607037},
															    {11799053, 12116660}};
//	public static final String[] CENTROMERE_BOUNDARIES = {"chr1:121,100,000-128,000,000", "chr2:91,000,000-95,700,000", "chr3:89,400,000-93,200,000", "chr4:48,700,000-52,400,000", "chr5:45,800,000-50,500,000", "chr6:58,400,000-63,400,000", "chr7:57,400,000-61,100,000", "chr8:43,200,000-48,100,000", "chr9:46,700,000-60,300,000", "chr10:38,800,000-42,100,000", "chr11:51,400,000-56,400,000", "chr12:33,200,000-36,500,000", "chr13:13,500,000-18,400,000", "chr14:13,600,000-19,100,000", "chr15:14,100,000-18,400,000", "chr16:34,400,000-40,700,000", "chr17:22,200,000-23,200,000", "chr18:15,400,000-17,300,000", "chr19:26,700,000-30,200,000", "chr20:25,700,000-28,400,000", "chr21:10,000,000-13,200,000", "chr22:9,600,000-16,300,000", "chrX:56,600,000-65,000,000", "chrY:11,300,001-12,500,000"};

	public static int[] parseUCSClocation(String str) {
    	int chr, start, stop;
    
    	// if (!str.startsWith("chr") || !str.contains(":") ||
    	// !str.contains("-")) {
    	if (!str.startsWith("chr")) {
    		System.err.println("Error - '"+str+"' is not a proper UCSC position");
    	}
    
    	try {
    		str = ext.replaceAllWith(str, ",", "");
    		if (str.contains(":-") || str.contains("--")) {
    			System.err.println("Warning - UCSC position '"+str+"' contains a negative position; returning whole chromosome");
    			chr = Integer.parseInt(str.substring(3, str.indexOf(":")));
    			start = -1;
    			stop = -1;
    		} else if (str.contains(":")) {
    			chr = Positions.chromosomeNumber(str.substring(3, str.indexOf(":")));
    			if (str.contains("-")) {
    				start = Integer.parseInt(str.substring(str.indexOf(":")+1, str.indexOf("-")));
    				stop = Integer.parseInt(str.substring(str.indexOf("-")+1));
    			} else {
    				start = stop = Integer.parseInt(str.substring(str.indexOf(":")+1));
    			}
    		} else if (str.endsWith("p")) {
    			chr = Integer.parseInt(str.substring(3, str.length()-1));
    			start = 0;
    			stop = CENTROMERE_MIDPOINTS[chr];
    		} else if (str.endsWith("q")) {
    			chr = Integer.parseInt(str.substring(3, str.length()-1));
    			start = CENTROMERE_MIDPOINTS[chr];
    			stop = CHROMOSOME_LENGTHS[chr];
    		} else {
    			chr = Integer.parseInt(str.substring(3));
    			start = -1;
    			stop = -1;
    		}
    	} catch (Exception e) {
    		System.err.println("Error - '"+str+"' is not a proper UCSC position");
    		return null;
    	}
    
    	return new int[] {chr, start, stop};
    }

	public static String getUCSCformat(String[] chr_pos) {
		int[] pos = new int[chr_pos.length];
		
		pos[0] = chromosomeNumber(chr_pos[0]);
		for (int i = 1; i < chr_pos.length; i++) {
			pos[i] = Integer.parseInt(chr_pos[i]);
		}
		
		return getUCSCformat(pos);
	}
	
	public static String getUCSCformat(int[] pos) {
    	if (pos.length!=2&&pos.length!=3) {
    		System.err.println("Error - not a valid UCSC position");
    		return null;
    	}
    	return "chr"+pos[0]+":"+pos[1]+(pos.length==3?"-"+pos[2]:"");
    }

	public static String getUCSClink(int[] pos) {
    	if (pos.length!=3) {
    		System.err.println("Error - not a valid UCSC position");
    		return null;
    	}
    	return "http://genome.ucsc.edu/cgi-bin/hgTracks?position="+getUCSCformat(pos);
    }

	public static String getUCSClinkInExcel(int[] pos) {
    	if (pos.length!=3) {
    		System.err.println("Error - not a valid UCSC position");
    		return null;
    	}
    	return "=HYPERLINK(\""+getUCSClink(pos)+"\", \""+getUCSCformat(pos)+"\")";
    }

	public static byte chromosomeNumber(String chromosome) {
    	byte chr = -1;
    	
    	if (chromosome.startsWith("chr")) {
    		chromosome = chromosome.substring(3);
    	}
    
    	if (chromosome.equals("XY") || chromosome.equals("PAR")) {
    		chr = 25;
    	} else if (chromosome.equalsIgnoreCase("X")) {
    		chr = 23;
    	} else if (chromosome.equalsIgnoreCase("Y")) {
    		chr = 24;
    	} else if (chromosome.equalsIgnoreCase("MT")||chromosome.equalsIgnoreCase("M")||chromosome.equalsIgnoreCase("Mito")) {
    		chr = 26;
    	} else if (chromosome.equalsIgnoreCase("un")||chromosome.equalsIgnoreCase("multi")||chromosome.equalsIgnoreCase("altonly")||chromosome.equalsIgnoreCase("noton")) {
    		chr = 0;
    	} else {
    		try {
    			chr = Byte.parseByte(chromosome);
    			if (chr<0||chr>27) {
    				System.err.println("Error - chromosome number '"+chromosome+"' is unassseptable");
    			}
    		} catch (NumberFormatException nfe) {
    			System.err.println("Error - '"+chromosome+"' is an invalid chromosome");
    		}
    	}
    
    	return chr;
    }

	public static String chromosomeNumberInverse(int chr) {
    	return Positions.CHR_CODES[chr];
    }

}
