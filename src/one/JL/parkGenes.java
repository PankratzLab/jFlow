package one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import common.Array;
import common.Logger;
import filesys.CNVariant;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.LocusSet;
import filesys.Segment;

public class parkGenes {

	
	public static void main(String[] args) {
		int numArgs = args.length;
		GeneTrack geneTrack = GeneTrack.load("N:/statgen/NCBI/RefSeq_hg19.gtrack", false);
		
		LocusSet<CNVariant> set = CNVariant.loadLocSet("D:/data/NGRC/cnvs/decentCalls_centromeresBroken.cnv", new Logger());
		ArrayList<CNVariant> found = new ArrayList<CNVariant>();
		ArrayList<String> genes = new ArrayList<String>();
		genes.add("HTRA2");
		genes.add("SNCA");
		genes.add("PINK1");
		genes.add("LRRK2");
		genes.add("ATP13A2");
		genes.add("PARK7");
		genes.add("PARK2");
		ArrayList<GeneData> genLocs = new ArrayList<GeneData>();
		for (String gene : genes) {
			if (geneTrack.lookupAllGeneData(gene).length == 0) {
				throw new IllegalArgumentException(gene);
			}
			GeneData[] d = geneTrack.lookupAllGeneData(gene);
			for (int i = 0; i < d.length; i++) {
				genLocs.add(d[i]);
			}
		}
		Segment[] segs = Segment.sortSegments(genLocs.toArray(new GeneData[genLocs.size()]));
		for (int i = 0; i < set.getLoci().length; i++) {
			if (Segment.overlapsAny(set.getLoci()[i], segs)) {
				found.add(set.getLoci()[i]);
			}
		}
		System.out.println(found.size());

		String out = "D:/data/NGRC/cnvs/pdGenes.cnv";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\tGENE");
			for (CNVariant cnv : found) {
				for (GeneData gene : genLocs) {
					if (cnv.overlaps(gene)) {
						writer.println(cnv.toPlinkFormat() + "\t" + gene.getGeneName());
						System.out.println("HDF");
					}
				}
			}
			writer.close();
		} catch (Exception e) {

		}
	}
}
