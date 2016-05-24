package kaput;

import java.io.*;

public class dumpibd {
	public dumpibd() throws IOException {
		PrintWriter writer = null;
		String chrom;
		writer = new PrintWriter(new FileWriter("batch.dumpibd"));

		writer.println("#/bin/sh\n\n");

		writer.println("cp struct111-.dat struct.dat");
		writer.println("cp ../../../chromo* .");
		writer.println("threeB");
		writer.println();

		for (int j = 1; j<=23; j++) {
			if (j<10) {
				chrom = "0"+j;
			} else {
				chrom = ""+j;
			}
			writer.println("java -classpath /home/npankrat/" + common.PSF.Java.GENVISIS + " park.bat.dat2loc map"+chrom+".dat");
			writer.println("echo -e \"pairs\\n3\\nload map"+chrom+".loc\\nprep re_chrom"+chrom+".pre\\nn\\nscan\\ndump ibd\\nmibd"+chrom+".dat\\nquit\\n\" | /software/bin/sibs");
		}

		writer.close();

		Runtime.getRuntime().exec("chmod +x batch.dumpibd");
	}

	public static void main(String[] args) throws IOException {
		try {
			new dumpibd();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
