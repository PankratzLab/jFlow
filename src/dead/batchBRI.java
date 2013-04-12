package dead;

import java.io.*;

public class batchBRI {
	// String[] SNPS = {"rs6714365", "rs2289235", "rs10933359", "rs934820"};
	// String[] SNP_POS = {"(-12.9)", "(+0)", "(+2.6)", "(+9.3)"};
	// String[] SNPS = {"rs10200894"};
	// String[] SNP_POS = {""};
	
	public static final String[] SNPS = {"haplotype-111X", "haplotype-222X", "haplotype-X1X2", "haplotype-X2X1"};

	public static final String[] SNP_POS = {"", "", "", ""};

	public static final String[] CUT3 = {"everyone", "no_parkin", "no_parkin_or_LRRK2"};

	public batchBRI(String jarPath, String genoPath) throws IOException {
		PrintWriter writer = null;
		// String opReps = "";
		// String opReps = " reps=50000";
		// String opReps = " reps=10000";
		String opReps = " reps=1000";

		// do counts on Windows
		writer = new PrintWriter(new FileWriter("spreadem.bat"));
		writer.println("cd BRI3");

		for (int i = 0; i<SNPS.length; i++) {
			writer.println("cd "+SNPS[i]+SNP_POS[i]);
			writer.println("cd probands");
			writer.println("copy probands_struct.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
			writer.println("move BRI_summary.out probands-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("cd ..");
			writer.println("cd plates_1-10");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("copy struct111-.dat struct.dat");
				writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
				writer.println("move BRI_summary.out 111-BRI_summary.out");
				writer.println("del struct.dat");
				writer.println("copy struct100-.dat struct.dat");
				writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
				writer.println("move BRI_summary.out 100-BRI_summary.out");
				writer.println("del struct.dat");
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd originalSFH");
			writer.println("copy 64SFH_struct.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
			writer.println("move BRI_summary.out 64SFH-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("cd ..");
			writer.println("cd plates_1-14");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("copy struct111-.dat struct.dat");
				writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
				writer.println("move BRI_summary.out 111-BRI_summary.out");
				writer.println("del struct.dat");
				writer.println("copy struct100-.dat struct.dat");
				writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
				writer.println("move BRI_summary.out 100-BRI_summary.out");
				writer.println("del struct.dat");
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd expandedSFH");
			writer.println("copy struct100-.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-cases.dat"+opReps);
			writer.println("move BRI_summary.out 100-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("cd ..");

			writer.println("cd controls");
			writer.println("copy CARES1-struct.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-controls.dat");
			writer.println("move BRI_summary.out CARES1-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("copy CARES2-struct.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-controls.dat");
			writer.println("move BRI_summary.out CARES2-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("copy LOAD-struct.dat struct.dat");
			writer.println("java -cp "+jarPath+" park.countBRI file="+genoPath+"BRI3-"+SNPS[i]+"-controls.dat");
			writer.println("move BRI_summary.out LOAD-BRI_summary.out");
			writer.println("del struct.dat");
			writer.println("cd ..");

			writer.println("cd ..");
		}
		writer.println("cd ..");
		writer.close();
		genoPath = "/work/npankrat/BRI3/";

		// do gist on unix
		writer = new PrintWriter(new FileWriter("gistem"));
		writer.println("cd BRI3");

		for (int i = 0; i<SNPS.length; i++) {
			writer.println("cd "+SNPS[i]+SNP_POS[i].replaceAll("\\(", "\\\\(").replaceAll("\\)", "\\\\)"));
			writer.println("cd plates_1-10");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("cp struct111-.dat struct.dat");
				writer.println(gistThis("111-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
				writer.println("cp struct100-.dat struct.dat");
				writer.println(gistThis("100-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd originalSFH");
			writer.println("cp 64SFH_struct.dat struct.dat");
			writer.println(gistThis("64SFH-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
			writer.println("cd ..");
			writer.println("cd plates_1-14");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("cp struct111-.dat struct.dat");
				writer.println(gistThis("111-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
				writer.println("cp struct100-.dat struct.dat");
				writer.println(gistThis("100-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd expandedSFH");
			writer.println("cp struct100-.dat struct.dat");
			writer.println(gistThis("100-", genoPath+"BRI3-"+SNPS[i]+"-cases.dat"));
			writer.println("cd ..");

			writer.println("cd ..");
		}
		writer.println("cd ..");
		writer.close();

		// run gist for all datasets on Windows
		writer = new PrintWriter(new FileWriter("gistAll.bat"));
		writer.println("cd BRI3");

		for (int i = 0; i<SNPS.length; i++) {
			writer.println("cd "+SNPS[i]+SNP_POS[i]);
			writer.println("cd plates_1-10");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("gist.exe 111-gist-2@233-1.dat > 111-gist-1.out");
				writer.println("gist.exe 111-gist-2@233-2.dat > 111-gist-2.out");
				writer.println("gist.exe 100-gist-2@233-1.dat > 100-gist-1.out");
				writer.println("gist.exe 100-gist-2@233-2.dat > 100-gist-2.out");
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd originalSFH");
			writer.println("gist.exe 64SFH-gist-2@233-1.dat > 64SFH-gist-1.out");
			writer.println("gist.exe 64SFH-gist-2@233-2.dat > 64SFH-gist-2.out");
			writer.println("cd ..");
			writer.println("cd plates_1-14");
			for (int j = 0; j<CUT3.length; j++) {
				writer.println("cd "+CUT3[j]);
				writer.println("gist.exe 111-gist-2@233-1.dat > 111-gist-1.out");
				writer.println("gist.exe 111-gist-2@233-2.dat > 111-gist-2.out");
				writer.println("gist.exe 100-gist-2@233-1.dat > 100-gist-1.out");
				writer.println("gist.exe 100-gist-2@233-2.dat > 100-gist-2.out");
				writer.println("cd ..");
			}
			writer.println("cd ..");
			writer.println("cd expandedSFH");
			writer.println("gist.exe 100-gist-2@233-1.dat > 100-gist-1.out");
			writer.println("gist.exe 100-gist-2@233-2.dat > 100-gist-2.out");
			writer.println("cd ..");

			writer.println("cd ..");
		}

		writer.println("cd ..");
		writer.close();

	}

	public String gistThis(String prefix, String filename) {
		return "twoB 2\n"+"jcp batch 8 2\n"+"./batch\n"+"jcp gist mut="+filename+" allele=1\n"+"jcp gist mut="+filename+" allele=2\n"+"cp "+filename.substring(0, filename.lastIndexOf("/")+1)+"gist.exe .\n"+"mv gist-2@233-1.dat "+prefix+"gist-2@233-1.dat\n"+"mv gist-2@233-2.dat "+prefix+"gist-2@233-2.dat\n"+"mv re_chrom02.pre "+prefix+"re_chrom02.pre\n"+"rm struct.dat chrom02.pre logfile\\ of\\ errors.out batch useful.opt allegro.log chrom02.lin.out chromf02.lin.out";
	}

	public static void main(String[] args) throws IOException {
		String jarpath = "c:\\park.jar";
		String genopath = "c:\\";

		String usage = "\n"+"park.genoOnIBD requires 3 arguments:\n"+"   (1) location of the park.jar file (i.e. jar="+jarpath+" (default), no spaces please)\n"+"   (2) path to the BRI-rs###-cases.dat files (i.e. genos="+genopath+" (default))\n"+"";
		int numArgs = args.length;

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("jar=")) {
				jarpath = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("genos=")) {
				genopath = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		System.out.println("Making a batch file utilizing "+jarpath+" and the files in "+genopath);
		try {
			new batchBRI(jarpath, genopath);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
