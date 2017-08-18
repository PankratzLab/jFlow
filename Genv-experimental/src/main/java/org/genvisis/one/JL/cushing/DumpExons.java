package org.genvisis.one.JL.cushing;

import java.util.HashSet;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.Segment;

public class DumpExons {

	public static void main(String[] args) {
		HashSet<String> genes = new HashSet<>();
		genes.add("CABLES1");
		genes.add("PRKAR1A");

		GeneTrack geneTrack = GeneTrack
				.load(Resources.genome(GENOME_BUILD.valueOf("HG19"), new Logger()).getGTrack().get());

		String top = "browser position chr1:14621-14712\ntrack name=\"Target Regions\"";

		StringBuilder out = new StringBuilder(top);
		for (GeneData[] geneDatas : geneTrack.getGenes()) {
			for (GeneData geneData : geneDatas) {
				if (genes.contains(geneData.getGeneName())) {
					for (int j = 0; j < geneData.getExonBoundaries().length; j++) {
						Segment exon = new Segment(geneData.getChr(), geneData.getExonBoundaries()[j][0],
								geneData.getExonBoundaries()[j][1]);
						out.append("\n"+exon.getChromosomeUCSC() + "\t" + exon.getStart() + "\t" + exon.getStop()
								+ "\tEXON_INDEX" + (j + 1) + "|" + geneData.getGeneName() + "|" + exon.getUCSClocation()
								);

					}
				}
			}
		}
		Files.write(out.toString(), "/Volumes/Beta/data/Cushings/QC/" + ArrayUtils.toStr(genes, "_") + ".exons.bed");
	}
}
