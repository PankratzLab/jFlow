package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.common.Files;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

/**
 * Helper utility for taking a CNV file and creating a regions list covering the samples of interest
 * and their parents (if present).
 */
public class CNVToFamilyRegions {

	public static void main(String... args) {
		final String cnvs = "cnvs";
		final String pedigree = "ped";

		CLI c = new CLI(PennCNVFamilies.class);
		c.addArg(cnvs, "Plink .cnv file from which to generate regions", true);
		c.addArg(pedigree, "Pedigree defining familial relationships for", true);

		c.parseWithExit(args);

		buildRegions(c.get(cnvs), c.get(pedigree));
	}

	private static void buildRegions(String cnvs, String pedigree) {
		final Pedigree ped = new Pedigree(null, pedigree, false);

		try {
			final String out = "familyRegions.txt";
			BufferedReader reader = Files.getAppropriateReader(cnvs);
			PrintWriter writer = Files.getAppropriateWriter(out);

			int[] idx = ext.indexFactors(new String[] {"FID", "IID", "CHR", "BP1", "BP2"},
			                             CNVariant.PLINK_CNV_HEADER, true, true);

			// Skip header
			reader.readLine();

			while (reader.ready()) {
				String[] line = reader.readLine().split("\t");
				String cFid = line[idx[0]];
				String cIid = line[idx[1]];
				String chr = line[idx[2]];
				String bp1 = line[idx[3]];
				String bp2 = line[idx[4]];

				int indIndex = ped.getIndIndex(cFid, cIid);
				int faIndex = ped.getIndexOfFaInIDs(indIndex);
				int moIndex = ped.getIndexOfMoInIDs(indIndex);

				String cDNA = ped.getiDNA(indIndex);
				StringBuilder comment = new StringBuilder(cDNA).append(", has father? ")
				                                               .append(faIndex > 0 ? "yes" : "no")
				                                               .append(" - has mother? ")
				                                               .append(moIndex > 0 ? "yes" : "no");

				Integer p1 = Integer.parseInt(bp1);
				Integer p2 = Integer.parseInt(bp2);
				int length = p2 - p1;
				p1 -= length;
				p2 += length;

				String pos = Positions.getUCSCformat(new String[] {chr, p1.toString(), p2.toString()});

				write(writer, cDNA, pos, comment.toString());

				if (faIndex > 0) {
					String faDNA = ped.getiDNA(faIndex);
					write(writer, faDNA, pos, "Father of " + cDNA);
				}

				if (moIndex > 0) {
					String moDNA = ped.getiDNA(moIndex);
					write(writer, moDNA, pos, "Mother of " + cDNA);
				}
			}

			System.out.println("Created regions file: " + out);
			writer.close();
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void write(PrintWriter writer, String sample, String pos, String comment) {
		writer.print(sample);
		writer.print("\t");
		writer.print(pos);
		writer.print("\t");
		writer.println(comment);
	}
}
