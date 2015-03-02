package gwas;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;

import bioinformatics.Sequence;
import common.Aliases;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.ext;

public class GeneScorePipeline {
	
	private static float DEFAULT_INDEX_THRESHOLD = (float)0.00000005;
	private static int DEFAULT_WINDOW_MIN_SIZE_PER_SIDE = 500000;// 500kb each side is technically a 1M window until the next hit region, but we now take this into consideration in the main algorithm
	private static float DEFAULT_WINDOW_EXTENSION_THRESHOLD = (float)0.00000005; // (float)0.00001;
	private static String[] DEFAULT_ADDL_ANNOT_VAR_NAMES = new String[0];
	
	private static final String[][] LINKERS = {
		{"MarkerName", "Marker", "SNP", "Variant", "VariantName"},
		Aliases.ALLELES[0],
		{"beta", "beta_SNP_add", "Effect", "b"} // aka Aliases.EFFECTS [which doesn't include 'b']
	};
	
	private String bfilePrefix;
	private String dataFile;
	
	private int bimChrIndex = 0;
	private int bimMkrIndex = 1;
	private int bimPosIndex = 3;
	private int bimA1Index = 4;
	private int bimA2Index = 5;
	
	private int dataMkrIndex = 0;
	private int dataA1Index = 3;
	private int dataBetaIndex = 6;
	
	private int hitsMkrIndex = 1;
	
	
	String subsetBIMDataFile;
	String hitWindowFile;
	String hitMkrDataFile;
	
	/*
	 * FILE FORMATS:
	 * 	    <bimFile> : .bim file format [TODO generalize for all plink-related marker info file formats]
	 *      <dataFile> : 
	 *      	1] MarkerName
	 *      	2] Allele1	
	 *      	3] Allele2	
	 *      	4] Freq.Allele1.HapMapCEU	
	 *      	5] b	
	 *      	6] SE	
	 *      	7] p	
	 *      	8] N
	 * 
	 */
	
	public GeneScorePipeline(String bFilePrefix, String dataFile) {
		this.bfilePrefix = bFilePrefix;
		this.dataFile = dataFile;
	}
	
	public void setBIMIndices(int mkrInd, int chrInd, int posInd) {
		bimChrIndex = chrInd;
		bimMkrIndex = mkrInd;
		bimPosIndex = posInd;
	}
	
	public void setDataIndices(int mkrInd, int a1Ind, int betaInd) {
		dataMkrIndex = mkrInd != -1 ? mkrInd : 0;
		dataA1Index = a1Ind != -1 ? a1Ind : 3;
		dataBetaIndex = betaInd != -1 ? betaInd : 6;
	}
	
	public void findHeaderIndices() {
		String[] header = Files.getHeaderOfFile(dataFile, null);
		int[] linkKeyIndices = ext.indexFactors(LINKERS, header, false, true, false, null, false);
		dataMkrIndex = linkKeyIndices[0];
		dataA1Index = linkKeyIndices[1];
		dataBetaIndex = linkKeyIndices[2];
	}

	private GeneScorePipeline crossFilterMarkerData() throws IOException {
		subsetBIMDataFile = ext.rootOf(dataFile, false) + ".bimData.xln";
		System.out.println("Cross-filtering data and .BIM files [ --> '" + subsetBIMDataFile + "']");
		BufferedReader bimReader;
		BufferedReader dataReader;
		PrintWriter dataWriter;
		HashMap<String, String[]> mkrsBim;
		
		mkrsBim = new HashMap<String, String[]>();
		bimReader = Files.getAppropriateReader(bfilePrefix + ".bim");
		String line = bimReader.readLine();
		int cntAmbig = 0;
		do {
			String[] parts = line.split("[\\s]+");
			String mkr = parts[bimMkrIndex];
			String a1 = parts[bimA1Index];
			String a2 = parts[bimA2Index];
			if (Sequence.validAllele(a1) && Sequence.validAllele(a2) && !a1.equals(Sequence.flip(a2))) {
				mkrsBim.put(mkr, parts);
			} else {
				cntAmbig++;
			}
		} while ((line = bimReader.readLine()) != null);
		bimReader.close();
		System.out.println("Found " + cntAmbig + " ambiguous markers (will be excluded)");
		dataReader = Files.getAppropriateReader(dataFile);
		dataWriter = new PrintWriter(subsetBIMDataFile);
		
		String[] dataHdrs = dataReader.readLine().split("[\\s]+");
		dataWriter.println("MarkerName\tChr\tPosition\t" + Array.toStr(Array.subArray(dataHdrs, 1))); //Allele1\tAllele2\tFreq.Allele1.HapMapCEU\tb\tSE\tp\tN
		while((line = dataReader.readLine()) != null) {
			String[] parts = line.split("[\\s]+");
			if (mkrsBim.containsKey(parts[dataMkrIndex])) {
				String[] bimParts = mkrsBim.get(parts[dataMkrIndex]);
				dataWriter.print(parts[dataMkrIndex]);
				dataWriter.print("\t");
				dataWriter.print(bimParts[bimChrIndex]);
				dataWriter.print("\t");
				dataWriter.print(bimParts[bimPosIndex]);
				dataWriter.print("\t");
				dataWriter.println(Array.toStr(Array.subArray(parts, 1)));
			}
		}
		dataWriter.flush();
		dataWriter.close();
		dataReader.close();
		
		return this;
	}
	
	private GeneScorePipeline runHitWindows(float indexThreshold, int windowMinSizePerSide, float windowExtensionThreshold, String[] additionalAnnotationVariableNames) {
		hitWindowFile = ext.rootOf(subsetBIMDataFile, false) + ".hits.out";
		System.out.println("Running hit window analysis [ --> '" + hitWindowFile + "']");
		String[][] results = HitWindows.determine(subsetBIMDataFile, indexThreshold, windowMinSizePerSide, windowExtensionThreshold, additionalAnnotationVariableNames);
		System.out.println("Found " + results.length + " hit windows");
		Files.writeMatrix(results, hitWindowFile, "\t");
		return this;
	}
	
	private GeneScorePipeline extractHitMarkerData() {
		hitMkrDataFile = ext.rootOf(subsetBIMDataFile, false) + ".subsetData.xln";
		System.out.println("Extracting data for hit window markers [ --> '" + hitMkrDataFile + "']");
		String[] hitMarkers = HashVec.loadFileToStringArray(hitWindowFile, true, new int[]{hitsMkrIndex}, false);
		HashSet<String> hitMrkSet = new HashSet<String>();
		for (String mkr : hitMarkers) {
			hitMrkSet.add(mkr);
		}
		
		int[] cols = new int[]{dataMkrIndex, dataA1Index, dataBetaIndex};
		String[][] bimData = HashVec.loadFileToStringMatrix(subsetBIMDataFile, true, cols, false);
		
		PrintWriter hitDataWriter = Files.getAppropriateWriter(hitMkrDataFile);
		hitDataWriter.println("MarkerName\tAllele1\tb");
		for (String[] markerData : bimData) {
			if (hitMrkSet.contains(markerData[dataMkrIndex])) {
				hitDataWriter.println(Array.toStr(markerData));
			}
		}
		hitDataWriter.flush();
		hitDataWriter.close();
		
		return this;
	}
	
	private void runPlink(boolean plink2) {
		System.out.print("Running plink command [ --> '");
		String dir = ext.parseDirectoryOfFile(hitMkrDataFile);
		String cmd = "plink" + (plink2 ? "2" : "") + " --noweb --bfile " + ext.rootOf(bfilePrefix+".bim", false) + " --score " + hitMkrDataFile;
		System.out.println(cmd + "']");
		/*boolean results = */CmdLine.run(cmd, dir);
	}
	
	private void writePlink(boolean plink2) {
		System.out.println("Writing plink command");
		String cmd = "plink" + (plink2 ? "2" : "") + " --noweb --bfile " + ext.rootOf(bfilePrefix+".bim", false) + " --score " + hitMkrDataFile;
		Files.write(cmd, ext.parseDirectoryOfFile(hitMkrDataFile) + "runPlink.sh");
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String bfile = "plink.bim";
		String dfile = "data.xln";
		float indexThreshold = DEFAULT_INDEX_THRESHOLD;
		int windowMinSizePerSide = DEFAULT_WINDOW_MIN_SIZE_PER_SIDE;
		float windowExtensionThreshold = DEFAULT_WINDOW_EXTENSION_THRESHOLD;
		
		int dataMkrIndex = -1;
		int dataA1Index = -1;
		int dataBetaIndex = -1;
		boolean findInHeader = false;
		boolean runPlink = false;
		boolean runPlink2 = false;
		
		String usage =  "\n" + 
						"lab.GeneScorePipeline requires 2+ arguments\n" + 
						"   (1) bfile; full path to plink .BIM file (i.e. bfile=" + bfile + " (default))\n" + 
						"   (2) dfile; full path to data file (i.e. dfile=" + dfile + " (default))\n" +
						"       OPTIONAL:\n" + 
						"   (3) p-value threshold for index SNPs (i.e. indexThresh=" + indexThreshold + " (default))\n" + 
						"   (4) minimum num bp per side of window (i.e. minWinSize=" + windowMinSizePerSide + " (default))\n" + 
						"   (5) p-value threshold to extend the window (i.e. winThresh=" + windowExtensionThreshold + " (default))\n" +
						"   (6) indices of necessary columns in data file [marker index, allele1 index, beta index] (i.e. cols=0,3,6 (default))\n" + 
						"   (7) search for indices of necessary columns in data file header (i.e. -findIndices (not the default))\n" + 
						"";
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bfile=")) {
				bfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dfile=")) {
				dfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("indexThresh=")) {
				indexThreshold = ext.parseFloatArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minWinSize=")) {
				windowMinSizePerSide = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("winThresh=")) {
				windowExtensionThreshold = ext.parseFloatArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("cols=")) {
				String[] colStrs = args[i].split("=")[1].split(",");
				if (colStrs.length < 3 || colStrs.length > 3) {
					System.err.println("Error - invalid argument: " + args[i]);
					System.err.println(usage);
					System.exit(1);
				}
				dataMkrIndex = Integer.parseInt(colStrs[0]);
				dataA1Index = Integer.parseInt(colStrs[1]);
				dataBetaIndex = Integer.parseInt(colStrs[2]);
				numArgs--;
			} else if (args[i].startsWith("-findIndices")) {
				findInHeader = true;
				numArgs--;
			} else if (args[i].startsWith("-runPlink")) {
				runPlink = true;
				runPlink2 = false;
				numArgs--;
			} else if (args[i].startsWith("-runPlink2")) {
				runPlink = false;
				runPlink2 = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0 || args.length == 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			GeneScorePipeline gsp = new GeneScorePipeline(bfile, dfile);
			if (findInHeader) {
				gsp.findHeaderIndices();
			} else if (dataMkrIndex != -1 || dataA1Index != -1 || dataBetaIndex != -1) {
				gsp.setDataIndices(dataMkrIndex, dataA1Index, dataBetaIndex);
			}
			gsp.crossFilterMarkerData();
			gsp.runHitWindows(indexThreshold, windowMinSizePerSide, windowExtensionThreshold, DEFAULT_ADDL_ANNOT_VAR_NAMES);
			gsp.extractHitMarkerData();
			if (runPlink || runPlink2) {
				gsp.runPlink(runPlink2);
			} else {
				gsp.writePlink(runPlink2);
			}
			System.out.println("Complete!");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}

