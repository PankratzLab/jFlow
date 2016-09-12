package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariant.CNVBuilder;

import org.genvisis.filesys.LocusSet;

public class CNV1000gConcordance {

	public static void main(String[] args) {
		String pfile = "/Volumes/Beta/data/1000G/1000GFormat.cnv";
		String dupFIle = "/Volumes/Beta/data/1000G/duplicates.cnv";
		String genFile = "/Volumes/Beta/data/1000G/genvisis.old.cnv";
		String finalCNVFile = "/Volumes/Beta/data/1000G/merge.cnv";

		Logger log = new Logger();
		if (!Files.exists(pfile)) {
			HashSet<String> samps = new HashSet<String>();

			samps.addAll(HashVec.loadFileToHashSet(dupFIle, new int[] { 0 }, "", false));

			samps.addAll(HashVec.loadFileToHashSet(dupFIle, new int[] { 1 }, "", false));
			log.reportTimeInfo(samps.size() + "");
			MarkerSet markerSet = MarkerSet.load("/Volumes/Beta/data/1000G/markers.ser", false);
			LocusSet<CNVariant> set = CNVariant
					.loadLocSet("/Volumes/Beta/data/1000G/GRCh37_hg19_variants_2015-07-23.txt.cnv", new Logger());
			int[][] indicesByChr = markerSet.getIndicesByChr();
			int num = 0;
			log.reportTimeInfo("HFD");
			PrintWriter writer = Files.getAppropriateWriter(pfile);
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < set.getLoci().length; i++) {
				// If the sample is in the genvisis project and is autosomal
				if (samps.contains(set.getLoci()[i].getIndividualID()) && set.getLoci()[i].getChr() < 23) {

					String[] markers = markerSet.getMarkersIn(set.getLoci()[i], indicesByChr);
					// If the number of markers >= 3
					if (markers != null && markers.length >= 10) {
						CNVBuilder builder = new CNVBuilder(set.getLoci()[i]);
						// Set the number of markers
						builder.numMarkers(markers.length);
						num++;
						writer.println(builder.build().toPlinkFormat());
					}
				}
				if (i % 10000 == 0) {
					log.reportTimeInfo(i + "\t" + num);
				}
			}
			writer.close();
		}
		LocusSet<CNVariant> set = CNVariant.loadLocSet(pfile, new Logger());

		HashSet<String> rep = CNVariant.getUniqueInds(set, log);
		HashSet<String> inds = new HashSet<String>();
		HashSet<String> finals = new HashSet<String>();

		for (String ind : rep) {
			inds.add(ind.split("\t")[0]);
		}
		ArrayList<String> comps = new ArrayList<String>();
		ArrayList<String> newSampDat = new ArrayList<String>();
		newSampDat.add("DNA\tFID\tIID");
		try {
			BufferedReader reader = Files.getAppropriateReader(dupFIle);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if (inds.contains(line[0]) || inds.contains(line[1])) {
					finals.add(line[0]);
					finals.add(line[1]);
					comps.add(Array.toStr(line));
					newSampDat.add(line[0] + "\t" + line[0] + "\t" + line[0]);
					newSampDat.add(line[1] + "\t" + line[1] + "\t" + line[1]);

				}
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Files.writeArrayList(comps, "/Volumes/Beta/data/1000G/finalDupSet.txt");
		Files.writeArrayList(newSampDat, "/Volumes/Beta/data/1000G/JLSampleData.txt");

		LocusSet<CNVariant> genset = CNVariant.loadLocSet(genFile, new Logger());
		PrintWriter writer = Files.getAppropriateWriter(finalCNVFile);
		writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
		for (int j = 0; j < set.getLoci().length; j++) {
			CNVBuilder builder = new CNVBuilder(set.getLoci()[j]);
			// builder.cn(1);
			writer.println(builder.build().toPlinkFormat());
		}
		for (int j = 0; j < genset.getLoci().length; j++) {
			if (finals.contains(genset.getLoci()[j].getIndividualID())) {
				CNVBuilder builder = new CNVBuilder(genset.getLoci()[j]);
				if (builder.getCn() == 4) {
					builder.cn(3);
				}
				if (builder.getCn() == 0) {
					builder.cn(1);
				}
				writer.println(builder.build().toPlinkFormat());
			}
		}
		writer.close();

		Hashtable<String, String> match = new Hashtable<String, String>();
		String[] one = HashVec.loadFileToStringArray(dupFIle, false, new int[] { 0 }, false);

	}

}
