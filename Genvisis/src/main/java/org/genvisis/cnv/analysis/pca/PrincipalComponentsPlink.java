package org.genvisis.cnv.analysis.pca;
//package cnv.analysis.pca;
//
//import common.Array;
//import common.Files;
//import common.Logger;
//import common.PSF;
//import common.ext;
//
///**
// * Class for dealing with Principle components in plink style Format, note that plink mds files contain an unused column prior to principal components
// *
// */
//public class PrincipalComponentsPlink extends PrincipalComponentsResiduals {
//
//	private PlinkSet plinkSet;
//	private PCLoading[] pcLoadings;
//	private String[] vars;
//
//	public PrincipalComponentsPlink(PlinkSet plinkSet, int numComponents, Logger log) {
//		super(plinkSet.getMds(), numComponents, log);
//
//		String[] samps = getSamplesToReport();
//		double[] pc1 = getBasisAt(2);
//		for (int i = 0; i < pc1.length; i++) {
//			System.out.println(samps[i] + "\t" + pc1[i]);
//		}
//	}
//
//	private static class PCLoading {
//		private String var;
//		private double[] data;
//
//		public PCLoading(String var, int numPCs) {
//			this.var = var;
//			this.data = new double[numPCs];
//		}
//
//	}
//
//	private static class PlinkSet {
//		public static String mdsExt = ".mds.mds";
//		private String mds,bim,bed,fam;
//		private Logger log;
//		private boolean validSet;
//
//		public PlinkSet(String mds, Logger log) {
//			super();
//			this.mds = mds;
//			this.log = log;
//			this.validSet = true;
//		}
//
//		public void init() {
//			if (!mds.endsWith(mdsExt)) {
//				log.reportTimeError("Expecting mds file " + mds + " to end with " + mdsExt);
//			} else {
//				String root = ext.rootOf(ext.rootOf(mdsExt, false), false);
//				String[] set = PSF.Plink.getPlinkBedBimFam(root);
//				this.bed = set[0];
//				this.bim = set[1];
//				this.fam = set[2];
//
//				if (Files.exists("", set) && Files.exists(mds)) {
//					log.reportTimeInfo("Found a complete plink set to work with");
//					validSet = false;
//				} else {
//					log.reportTimeError("Missing one of the following files " + Array.toStr(set, "\n") + "\n" + mds);
//				}
//			}
//		}
//		
//		
//
//		public String getMds() {
//			return mds;
//		}
//
//	}
//
//	public static void test(String fullPathToPlinkRoot, int numComponents, Logger log) {
//		PlinkSet plinkSet = new PlinkSet(fullPathToPlinkRoot + PlinkSet.mdsExt, log);
//		PrincipalComponentsPlink pcPlink = new PrincipalComponentsPlink(plinkSet, numComponents, log);
//	}
//
//	public static void main(String[] args) {
//		int numArgs = args.length;
//		String fullPathToPlinkRoot = "/home/usr/plink/plink";
//		int numComponents = 3;
//		String logfile = null;
//		Logger log;
//
//		String usage = "\n" + "cnv.analysis.pca.PrincipalComponentsPlink requires 0-1 arguments\n";
//		usage += "   (1) full path to plink root (i.e. plinkRoot=" + fullPathToPlinkRoot + " (default))\n" + "";
//		usage += "   (2) number of components  (i.e. numComponents=" + numComponents + " (default))\n" + "";
//
//		for (int i = 0; i < args.length; i++) {
//			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
//				System.err.println(usage);
//				System.exit(1);
//			} else if (args[i].startsWith("plinkRoot=")) {
//				fullPathToPlinkRoot = args[i].split("=")[1];
//				numArgs--;
//			} else if (args[i].startsWith("numComponents=")) {
//				numComponents = ext.parseIntArg(args[i]);
//				numArgs--;
//			} else if (args[i].startsWith("log=")) {
//				logfile = args[i].split("=")[1];
//				numArgs--;
//			} else {
//				System.err.println("Error - invalid argument: " + args[i]);
//			}
//		}
//		if (numArgs != 0) {
//			System.err.println(usage);
//			System.exit(1);
//		}
//		try {
//			log = new Logger(fullPathToPlinkRoot + ".genvisis.log");
//			test(fullPathToPlinkRoot, numComponents + 1, log);
//			;
//			// parse(filename, log);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//	}
//
//}
