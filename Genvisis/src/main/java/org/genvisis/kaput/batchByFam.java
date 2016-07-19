package kaput;

import java.io.*;

public class batchByFam {

	public static boolean multiples = true;

	// public static int[] dist = {0, 222, 194, 187, 158, 150, 146, 147, 136,
	// 132, 146, 127, 137, 108, 108, 85, 111, 109, 114, 114, 105, 105, 68, 79};
	public static int[] dist = {0, 274, 261, 225, 213, 195, 187, 182, 164, 161, 178, 147, 166, 111, 125, 123, 131, 125, 120, 102, 99, 57, 46, 1};

	@SuppressWarnings("resource")
	public batchByFam() throws IOException {
		PrintWriter writer = null, writer1 = null, writer2 = null, writer3 = null, writer4 = null, writer5 = null, writer6 = null;
		int chromosome = 0;

		if (multiples) {
			writer1 = new PrintWriter(new FileWriter("batch.phenos.1"));
			writer2 = new PrintWriter(new FileWriter("batch.phenos.2"));
			writer3 = new PrintWriter(new FileWriter("batch.phenos.3"));
			writer4 = new PrintWriter(new FileWriter("batch.phenos.4"));
			writer5 = new PrintWriter(new FileWriter("batch.phenos.5"));
			writer6 = new PrintWriter(new FileWriter("batch.phenos.6"));
			writer1.println("#/bin/sh");
			writer1.println();
			writer1.println("sleep 30");
			writer1.println();
			writer2.println("#/bin/sh");
			writer2.println();
			writer2.println("sleep 30");
			writer2.println();
			writer3.println("#/bin/sh");
			writer3.println();
			writer3.println("sleep 30");
			writer3.println();
			writer4.println("#/bin/sh");
			writer4.println();
			writer4.println("sleep 30");
			writer4.println();
			writer5.println("#/bin/sh");
			writer5.println();
			writer5.println("sleep 30");
			writer5.println();
			writer6.println("#/bin/sh");
			writer6.println();
			writer6.println("sleep 30");
			writer6.println();
		} else {
			writer = new PrintWriter(new FileWriter("batch.phenos"));
			writer.println("#/bin/sh");
			writer.println();
			writer.println();
			writer.println("sleep 60");
			writer.println();
		}

		for (int i = 1; i<=330; i++) {
			if (multiples) {
				writer = null;
			}
			if (i>=1&&i<=50) {
				writer = writer1;
			}
			if (i>=51&&i<=100) {
				writer = writer2;
			}
			if (i>=101&&i<=150) {
				writer = writer3;
			}
			if (i>=151&&i<=210) {
				writer = writer4;
			}
			if (i>=211&&i<=270) {
				writer = writer5;
			}
			if (i>=271&&i<=330) {
				writer = writer6;
			}
			if (writer==null) {
				System.err.println("Error: Not all chromosomes were told which file to be in.");
			}

			for (int j = 0; j<2; j++) {
				switch (j) {
				case 0:
					chromosome = 15;
					break;
				case 1:
					chromosome = 21;
					break;
				default:
					System.err.println("Crap: Error with the whole which chromosome business");
					break;
				}

				// writer.println("cp solar_phen-age4549comp_mhgt.dat chrom"+
				// chrom);
				// writer.println("cp solar_phen-age4549comp_mbmi.dat chrom"+
				// chrom);
				// writer.println("cp solar_phen-age4549comp_mchol.dat chrom"+
				// chrom);
				// writer.println("cp solar_phen-age4549comp_msbp+comp_mbmi.dat
				// chrom"+ chrom);
				// writer.println("cp solar_phen-age4549comp_mwgt+comp_mhgt.dat
				// chrom"+ chrom);
				// writer.println("cp solar_phen-age5054comp_mgluc+comp_mbmi.dat
				// chrom"+ chrom);

				// writer.println("cp solar_phen-rate.dat chrom"+ chrom);
				// writer.println("cp solar_phen-pre_trt.dat chrom"+ chrom);
				// writer.println("cp solar_phen-post_trt.dat chrom"+ chrom);
				// writer.println("cp solar_phen-1stage.dat chrom"+ chrom);
				// writer.println("cp solar_phen-2ndage.dat chrom"+ chrom);
				writer.println("cp solar_phen-avgbp.dat fam"+i+"/chrom"+chromosome+"/");

				// writer.println("cp solar_phen-mdrink.dat chrom"+ chrom);
				// writer.println("cp solar_phen-mcpd.dat chrom"+ chrom);

				// writer.println("cp solar_phen-procam_auc.dat chrom"+ chrom);
				// writer.println("cp solar_phen-noage_procam_auc.dat chrom"+
				// chrom);
				// writer.println("cp solar_phen-pca1_auc.dat chrom"+ chrom);
				// writer.println("cp solar_phen-pca2_auc.dat chrom"+ chrom);

				writer.println("cd fam"+i+"/chrom"+chromosome+"/");
				writer.println("echo -e \""
				// +"load pedigree solar_ped."+chrom+"\\nload
				// freq solar_freq."+chrom+"\\nload marker
				// solar_marker."+chrom+"\\nibddir .\\nverbosity
				// min\\nibd\\nload map map."+chrom+"\\nibddir
				// .\\nmibddir .\\nmibd 0 "+dist[i]+"
				// 1\\nmibddir .\\n"

				// +"automodel solar_phen-rate.dat
				// rate\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-pre_trt.dat
				// pre_trt\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-post_trt.dat
				// post_trt\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-1stage.dat
				// 1stage\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-2ndage.dat
				// 2ndage\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				+"automodel solar_phen-avgbp.dat avgbp\\npolygenic -screen\\nmibddir .\\nchromosome "+chromosome+"\\ninterval 1\\nmultipoint -overwrite\\n"

				// +"automodel solar_phen-age4549comp_mhgt.dat
				// comp_mhgt\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-age4549comp_mbmi.dat
				// comp_mbmi\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-age4549comp_mchol.dat
				// comp_mchol\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel
				// solar_phen-age4549comp_msbp+comp_mbmi.dat
				// comp_msbp\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel
				// solar_phen-age4549comp_mwgt+comp_mhgt.dat
				// comp_mwgt\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel
				// solar_phen-age5054comp_mgluc+comp_mbmi.dat
				// comp_mgluc\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"

				// +"automodel solar_phen-mdrink.dat
				// mdrink\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-mcpd.dat
				// mcpd\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"

				// +"automodel solar_phen-procam_auc.dat
				// procam_auc\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-noage_procam_auc.dat
				// noage_procam_auc\\npolygenic
				// -screen\\nmibddir .\\nchromosome
				// "+i+"\\ninterval 1\\nmultipoint
				// -overwrite\\n"
				// +"automodel solar_phen-pca1_auc.dat
				// pca1_auc\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"
				// +"automodel solar_phen-pca2_auc.dat
				// pca2_auc\\npolygenic -screen\\nmibddir
				// .\\nchromosome "+i+"\\ninterval
				// 1\\nmultipoint -overwrite\\n"

				+"quit\\n\" | solar > fam"+i+".log");
				writer.println("cd ../..");
				// writer.println("lp -dib139 chrm"+ chrom +"-p4.ps");

			}
			writer.println();
		}

		if (multiples) {
			writer1.close();
			writer2.close();
			writer3.close();
			writer4.close();
			writer5.close();
			writer6.close();
		} else {
			writer.close();
		}

	}

	public static void main(String[] args) throws IOException {
		new batchByFam();
	}
}

// BufferedReader reader = null;
// PrintWriter writer = null;
// String temp, ID, chrome, blank="0";
// StringTokenizer st;
// Hashtable hash = new Hashtable();
//
// reader = new BufferedReader(new FileReader("dir.dir"));
// writer = new PrintWriter(new FileWriter("batch"));
//
// while (reader.ready()) {
// temp = reader.readLine();
// temp = temp.substring(0, temp.length() - 7);
// writer.println("gunzip -c "+temp+".tar.gz | tar xvf -");
//
// writer.println("cd "+temp);
// writer.println("rm -r chrom01");
// writer.println("rm -r chrom02");
// writer.println("rm -r chrom03");
// writer.println("rm -r chrom04");
// writer.println("rm -r chrom05");
// writer.println("rm -r chrom06");
// writer.println("rm -r chrom07");
// writer.println("rm -r chrom08");
// writer.println("rm -r chrom09");
// writer.println("rm -r chrom10");
// writer.println("rm -r chrom11");
// writer.println("rm -r chrom12");
// writer.println("rm -r chrom13");
// writer.println("rm -r chrom14");
// writer.println("rm -r chrom16");
// writer.println("rm -r chrom17");
// writer.println("rm -r chrom18");
// writer.println("rm -r chrom19");
// writer.println("rm -r chrom20");
// writer.println("rm -r chrom22");
// writer.println("cd ..");
//
// writer.println("tar -cvf "+temp+".tar "+temp);
// writer.println("rm "+temp+".tar.gz");
// writer.println("gzip "+temp+".tar");
// writer.println("rm -r "+temp);
// writer.println();
// }
// reader.close();
// writer.close();
