package org.genvisis.link.init;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import org.genvisis.common.ext;
import org.genvisis.link.Markers;

public class heterozygosity {

	public heterozygosity() throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		StringTokenizer st;
		String chrome;
		String[] markers, hetExp;
		double homos, homes;
		int numMarkers, numAlleles;

		writer = new PrintWriter(new FileWriter("heterozygosity summary.out"));
		writer.println("marker\tExp\tEmp");
		for (int i = 1; i <= 23; i++) {
			chrome = (i < 10) ? "0" + i : "" + i;
			try {
				// reader = new BufferedReader(new
				// FileReader("/home/npankrat/park/00masters/map"+chrome+".dat"));
				reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));

				st = new StringTokenizer(reader.readLine());
				numMarkers = Integer.valueOf(st.nextToken()).intValue() - 1;
				markers = new String[numMarkers];
				hetExp = new String[numMarkers];

				for (int j = 0; j < 6; j++) {
					reader.readLine();
				}
				if (i == 23) {
					reader.readLine();
				}
				for (int j = 0; j < numMarkers; j++) {
					st = new StringTokenizer(reader.readLine());
					st.nextToken();
					numAlleles = Integer.valueOf(st.nextToken()).intValue();
					st.nextToken();
					markers[j] = st.nextToken();
					st = new StringTokenizer(reader.readLine());
					homos = 0;
					for (int k = 0; k < numAlleles; k++) {
						homes = Double.valueOf(st.nextToken()).doubleValue();
						homos += homes * homes;
					}
					hetExp[j] = ext.formDeci(1 - homos, 2, true);
				}
				Markers.order(markers, true);
				reader = new BufferedReader(new FileReader("markerMap.dat"));
				reader.readLine();
				for (int j = 0; j < numMarkers; j++) {
					st = new StringTokenizer(reader.readLine());
					st.nextToken();
					writer.println(markers[j] + "\t" + hetExp[j] + "\t" + st.nextToken());
					if (j != numMarkers - 1) {
						reader.readLine();
					}
				}

			} catch (Exception e) {
				System.out.println("Problem with chromosome " + i);
			}

		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		try {
			new heterozygosity();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
