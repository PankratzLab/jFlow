//package cnv.hmm;
//
//import java.util.ArrayList;
//import java.util.Hashtable;
//
//import org.apache.commons.math3.random.MersenneTwister;
//
//import common.Array;
//import cnv.filesys.MarkerSet;
//import cnv.filesys.Project;
//import cnv.filesys.Sample;
//import cnv.hmm.CBSUtils.SegmentSplitUndo;
//import cnv.hmm.CBSUtils.SegmentationMethod;
//import filesys.Segment;
//
///**
// * @author lane0212 Genvisis version of Circular binary segmentation: might be applicable to both array and NGS...
// * 
// *
// */
//public class CircDuCnv {
//
//	private static final double DefaultAlpha = 0.01;
//	private static final String DEFAULT_DATA_TYPE = "logratio";
//	private String InputBinPath;
//	private String DataType;
//	private static final int idxChr = 0, idxStart = 1, idxEnd = 2;
//	private int idxScore = 3;
//	private Hashtable<String, int[]> StartByChr;
//	private Hashtable<String, int[]> EndByChr;
//	private Hashtable<String, double[]> ScoreByChr;
//	private GenomeSegmentationResults SegmentationResults;
//	public SegmentSplitUndo UndoMethod = SegmentSplitUndo.None;
//	private String ForbiddenIntervalBedPath = null;
//
//	// dataType: "logratio" (aCGH, ROMA, etc.) or "binary" (LOH)
//	public CircDuCnv(String inputBinPath, String dataType, int idxScore, Hashtable<String, int[]> startByChr, Hashtable<String, int[]> endByChr, Hashtable<String, double[]> scoreByChr, SegmentSplitUndo undoMethod, String forbiddenIntervalBedPath) {
//		super();
//		InputBinPath = inputBinPath;
//		DataType = dataType;
//		this.idxScore = idxScore;
//		StartByChr = startByChr;
//		EndByChr = endByChr;
//		ScoreByChr = scoreByChr;
//		UndoMethod = undoMethod;
//		ForbiddenIntervalBedPath = forbiddenIntervalBedPath;
//	}
//
//	public void SegmentGenome(String outPath, SegmentationMethod method) {
//		switch (method) {
//		case Wavelets:
//		default:// use Wavelets if CBS is not selected
//			System.out.println("Running Wavelet Partitioning");
//			break;
//		case CBS:
//			System.out.println("Running CBS Partitioning");
//			// CBS( 2);
//			break;
//		}
//		// this.WriteCanvasPartitionResults(outPath);
//	}
//
//	private GenomeSegmentationResults GetDummySegmentationResults() {
//		GenomeSegmentationResults results = new GenomeSegmentationResults(new Hashtable<String, Segment[]>());
//		return results;
//	}
//
//	// / <summary>
//	// / CBS: circular binary segmentation porting the R function segment in DNAcopy
//	// / </summary>
//	// / <param name="alpha">Now in this.Alpha</param>
//	// / <param name="nPerm"></param>
//	// / <param name="pMethod">"hybrid" or "perm"</param>
//	// / <param name="minWidth"></param>
//	// / <param name="kMax"></param>
//	// / <param name="nMin"></param>
//	// / <param name="eta"></param>
//	// / <param name="sbdry"></param>
//	// / <param name="trim"></param>
//	// / <param name="undoSplit">"none" or "prune" or "sdundo"; now in this.UndoMethod</param>
//	// / <param name="undoPrune"></param>
//	// / <param name="undoSD"></param>
//	// / <param name="verbose"></param>
//	// private void CBS(uint nPerm = 10000, string pMethod = "hybrid", int minWidth = 2, int kMax = 25,
//	// uint nMin = 200, double eta = 0.05, uint[] sbdry = null, double trim = 0.025,
//	// double undoPrune = 0.05, double undoSD = 3, int verbose = 1)
//	// {
//
//	private void CBS(int nPerm, String pMethod, int minWidth, int kMax, int nMin, double eta, int[] sbdry, double trim, double undoPrune, double undoSD, int verbose) {
//		if (minWidth < 2 || minWidth > 5) {
//			System.out.println("Minimum segment width should be between 2 and 5");
//			System.exit(1);
//		}
//		if (nMin < 4 * kMax) {
//			System.out.println("nMin should be >= 4 * kMax");
//			System.exit(1);
//		}
//		if (sbdry == null) {
//			sbdry =new int[1];
//			CBSUtils.ComputeBoundary(nPerm, DefaultAlpha, eta, sbdry);
//		}
//
//		Hashtable<String, int[]> inaByChr = new Hashtable<String, int[]>();
//		Hashtable<String, double[]> finiteScoresByChr = new Hashtable<String, double[]>();
//
//
//		for (String key : ScoreByChr.keySet()) {
//
//			int[] ina = CBSUtils.GetFiniteIndices(ScoreByChr.get(key));
//			System.out.println("ONfdsf chsdfsdfr " + key + "\t" + ina.length+"\t"+ScoreByChr.get(key).length);
//
//			double[] scores = null;
//			if (ina.length == ScoreByChr.get(key).length) {
//				scores = ScoreByChr.get(key);
//			} else {
//				scores = Array.subArray( ScoreByChr.get(key), ina);
//			}
//			finiteScoresByChr.put(key, scores);
//			inaByChr.put(key, ina);
//		}
//
//		// Quick sanity-check: If we don't have any segments, then return a dummy result.
//		int n = 0;
//		for (String list : finiteScoresByChr.keySet()) {
//			n += finiteScoresByChr.get(list).length;
//			System.out.println("ONf score f chr "+list+"\t"+n);
//
//		}
//		if (n == 0) {
//			this.SegmentationResults = this.GetDummySegmentationResults();
//			return;
//		}
//
//		double trimmedSD = Math.sqrt(CBSUtils.TrimmedVariance(finiteScoresByChr, trim));
//		Hashtable<String, Segment[]> segmentByChr = new Hashtable<String, Segment[]>();
//		MersenneTwister seedGenerator = new MersenneTwister(0);
//		// when parallelizing we need an RNG for each chromosome to get deterministic results
//		Hashtable<String, MersenneTwister> perChromosomeRandom = new Hashtable<String, MersenneTwister>();
//		for (String chr : ScoreByChr.keySet()) {
//			System.out.println("ON twister chr "+chr);
//
//			perChromosomeRandom.put(chr, new MersenneTwister(seedGenerator.nextInt()));
//		}
//
//		for (String chr : ScoreByChr.keySet()) {
//			System.out.println("ON chr "+chr);
//			int[] ina = inaByChr.get(chr);
//			int[] lengthSeg = new int[0];
//			double[] segmentMeans = new double[0];
//
//			lengthSeg = CBSUtils.ChangePoints(ScoreByChr.get(chr), sbdry, lengthSeg, segmentMeans, perChromosomeRandom.get(chr),
//
//			DataType, DefaultAlpha, nPerm, pMethod, minWidth, kMax, nMin, trimmedSD, UndoMethod, undoPrune, undoSD, verbose, 100, 1E-6);
//
//			SegmentCBS[] segments = new SegmentCBS[lengthSeg.length];
//			int cs1 = 0, cs2 = -1; // cumulative sum
//			for (int i = 0; i < lengthSeg.length; i++) {
//				System.out.println(lengthSeg.length+"\t"+lengthSeg[i]);
//				cs2 += lengthSeg[i];
//				int start = ina[cs1];
//				int end = ina[Math.min(cs2,ina.length-1)];
//				segments[i] = new SegmentCBS(chr, StartByChr.get(chr)[start], EndByChr.get(chr)[end], lengthSeg[i], 0);
//				System.out.println(segments[i].getUCSClocation());
//				cs1 += lengthSeg[i];
//			}
//
//		}
//
//		this.SegmentationResults = new GenomeSegmentationResults(segmentByChr);
//	}
//
//	private static class SegmentCBS extends Segment {
//		/**
//		 * 
//		 */
//		private static final long serialVersionUID = 1L;
//		private int nMarkers;
//		private double mean;
//
//		public SegmentCBS(String chr, int idxStart, int stop, int nMarkers, double mean) {
//			super(chr, idxStart, stop);
//			this.nMarkers = nMarkers;
//			this.mean = mean;
//		}
//	}
//
//	private class GenomeSegmentationResults {
//		private Hashtable<String, Segment[]> SegmentByChr;
//
//		public GenomeSegmentationResults(Hashtable<String, Segment[]> segmentByChr) {
//			this.SegmentByChr = segmentByChr;
//		}
//	}
//
//	private static void test(Project proj) {
//		Sample samp = proj.getFullSampleFromRandomAccessFile(proj.getSamples()[0]);
//		int nPerm = 10000;
//		String pMethod = "hybrid";
//		int minWidth = 2;
//		int kMax = 25;
//
//		int nMin = 200;
//		double eta = 0.05;
//		int[] sbdry = null;
//		double trim = 0.025;
//
//		double undoPrune = 0.05;
//		double undoSD = 3;
//		int verbose = 1;
//
//		Hashtable<String, int[]> StartByChr = new Hashtable<String, int[]>();
//		Hashtable<String, int[]> EndByChr = new Hashtable<String, int[]>();
//		Hashtable<String, double[]> ScoreByChr = new Hashtable<String, double[]>();
//		MarkerSet markerSet = proj.getMarkerSet();
//		double[] lrrs = Array.toDoubleArray(samp.getLRRs());
//		int[][] indices = proj.getMarkerSet().getIndicesByChr();
//
//		for (int i = 0; i < indices.length; i++) {
//			StartByChr.put(i + "", Array.subArray(markerSet.getPositions(), indices[i]));
//			EndByChr.put(i + "", Array.subArray(markerSet.getPositions(), indices[i]));
//			ScoreByChr.put(i + "", Array.subArray(lrrs, indices[i]));
//
//		}
//		CircDuCnv circDuCnv = new CircDuCnv(null, DEFAULT_DATA_TYPE, 3, StartByChr, EndByChr, ScoreByChr, SegmentSplitUndo.None, null);
//		circDuCnv.CBS(nPerm, pMethod, minWidth, kMax, nMin, eta, sbdry, trim, undoPrune, undoSD, verbose);
//	
//		// private void CBS(int nPerm, String pMethod, int minWidth, int kMax, int nMin, double eta, int[] sbdry, double trim, double undoPrune, double undoSD, int verbose) {
//
//	}
//public static void main(String[] args) {
//	int numArgs = args.length;
//	String filename = "CircDuCnv.dat";
//	String logfile = null;
//
//	String usage = "\n" + "cnv.hmm.CircDuCnv requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";
//
//	for (int i = 0; i < args.length; i++) {
//		if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
//			System.err.println(usage);
//			System.exit(1);
//		} else if (args[i].startsWith("file=")) {
//			filename = args[i].split("=")[1];
//			numArgs--;
//		} else if (args[i].startsWith("log=")) {
//			logfile = args[i].split("=")[1];
//			numArgs--;
//		} else {
//			System.err.println("Error - invalid argument: " + args[i]);
//		}
//	}
//	if (numArgs != 0) {
//		System.err.println(usage);
//		System.exit(1);
//	}
//	try {
//		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
//		test(proj);
//	} catch (Exception e) {
//		e.printStackTrace();
//	}
//}
//	// // / <summary>
//	// // / Assume that the rows are sorted by the start position and ascending order
//	// // / </summary>
//	// private void ReadBEDInput()
//	// {
//	// try
//	// {
//	// Dictionary<string, List<uint>> startByChr = new Dictionary<string, List<uint>>(),
//	// endByChr = new Dictionary<string, List<uint>>();
//	// Dictionary<string, List<double>> scoreByChr = new Dictionary<string, List<double>>();
//	// // Create an instance of StreamReader to read from a file.
//	// // The using statement also closes the StreamReader.
//	// using (GzipReader reader = new GzipReader(this.InputBinPath))
//	// {
//	// string line;
//	// string[] tokens;
//	// while ((line = reader.ReadLine()) != null)
//	// {
//	// tokens = line.Split('\t');
//	// string chr = tokens[Segmentation.idxChr].Trim();
//	// if (!startByChr.ContainsKey(chr))
//	// {
//	// startByChr.Add(chr, new List<uint>());
//	// endByChr.Add(chr, new List<uint>());
//	// scoreByChr.Add(chr, new List<double>());
//	// }
//	// startByChr[chr].Add(Convert.ToUInt32(tokens[Segmentation.idxStart].Trim()));
//	// endByChr[chr].Add(Convert.ToUInt32(tokens[Segmentation.idxEnd].Trim()));
//	// scoreByChr[chr].Add(Convert.ToDouble(tokens[this.idxScore].Trim()));
//	// }
//	// foreach (string chr in startByChr.Keys)
//	// {
//	// this.StartByChr[chr] = startByChr[chr].ToArray();
//	// this.EndByChr[chr] = endByChr[chr].ToArray();
//	// this.ScoreByChr[chr] = scoreByChr[chr].ToArray();
//	// }
//	//
//	// }
//	// }
//	// catch (Exception e)
//	// {
//	// Console.Error.WriteLine("File {0} could not be read:", this.InputBinPath);
//	// Console.Error.WriteLine(e.Message);
//	// Environment.Exit(1);
//	// }
//	// }
//	//
//	// private void WriteCanvasPartitionResults(string outPath)
//	// {
//	// Dictionary<string, bool> starts = new Dictionary<string, bool>();
//	// Dictionary<string, bool> stops = new Dictionary<string, bool>();
//	//
//	// foreach (string chr in SegmentationResults.SegmentByChr.Keys)
//	// {
//	// for (int segmentIndex = 0; segmentIndex < SegmentationResults.SegmentByChr[chr].Length; segmentIndex++)
//	// {
//	// Segment segment = SegmentationResults.SegmentByChr[chr][segmentIndex];
//	// starts[chr + ":" + segment.start] = true;
//	// stops[chr + ":" + segment.end] = true;
//	// }
//	// }
//	//
//	// Dictionary<string, List<GenomicBin>> ExcludedIntervals = new Dictionary<string, List<GenomicBin>>();
//	// if (!string.IsNullOrEmpty(ForbiddenIntervalBedPath))
//	// {
//	// ExcludedIntervals = CanvasCommon.Utilities.LoadBedFile(ForbiddenIntervalBedPath);
//	// }
//	//
//	// using (GzipWriter writer = new GzipWriter(outPath))
//	// {
//	// int segmentNum = -1;
//	//
//	// foreach (string chr in StartByChr.Keys)
//	// {
//	// List<GenomicBin> excludeIntervals = null;
//	// if (ExcludedIntervals.ContainsKey(chr)) excludeIntervals = ExcludedIntervals[chr];
//	// int excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
//	// uint previousBinEnd = 0;
//	// for (int pos = 0; pos < StartByChr[chr].Length; pos++)
//	// {
//	// uint start = StartByChr[chr][pos];
//	// uint end = EndByChr[chr][pos];
//	// bool newSegment = false;
//	// string key = chr + ":" + start;
//	// if (starts.ContainsKey(key))
//	// {
//	// newSegment = true;
//	// }
//	//
//	// if (excludeIntervals != null)
//	// {
//	// while (excludeIndex < excludeIntervals.Count && excludeIntervals[excludeIndex].Stop < previousBinEnd) excludeIndex++;
//	// if (excludeIndex < excludeIntervals.Count)
//	// {
//	// // Note: forbiddenZoneMid should never fall inside a bin, becuase these intervals were already excluded
//	// // from consideration during the call to CanvasBin.
//	// int forbiddenZoneMid = (excludeIntervals[excludeIndex].Start + excludeIntervals[excludeIndex].Stop) / 2;
//	// if (previousBinEnd < forbiddenZoneMid && end >= forbiddenZoneMid) newSegment = true;
//	// }
//	// }
//	// if (newSegment) segmentNum++;
//	// writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}", chr, start, end, ScoreByChr[chr][pos], segmentNum));
//	// previousBinEnd = end;
//	// }
//	// }
//	// }
//	// }
//}
