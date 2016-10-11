package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

public class simParseAllegro {
	public simParseAllegro(int start, int finish) throws IOException {
		BufferedReader input;
		PrintWriter writer = new PrintWriter(new FileWriter("summary.out"));
		String pos, maxPos;
		StringTokenizer st;
		float lod, maxLod;

		int repNum = start;
		try {
			while (repNum <= finish) {
				input = new BufferedReader(new FileReader("output-" + repNum + ".prn"));
				maxPos = "";
				maxLod = -1000;

				input.readLine();
				while (input.ready()) {
					st = new StringTokenizer(input.readLine());
					pos = st.nextToken();
					lod = Float.valueOf(st.nextToken()).floatValue();
					if (lod > maxLod) {
						maxLod = lod;
						maxPos = pos;
					}
				}

				writer.println(maxPos + "\t" + maxLod);
				repNum++;
				input.close();
			}
		} catch (Exception e) {
			System.err.println("Stopped prematurely at rep " + repNum + ".");
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			try {
				new simParseAllegro(1, 1000);
			} catch (Exception e) {
				e.printStackTrace();
			}
			// System.out.println("Expecting 2 arguments. The alpha and
			// omega.");
		} else {
			try {
				new simParseAllegro(Integer.valueOf(args[0]).intValue(),
														Integer.valueOf(args[1]).intValue());
				// new parseoutputs(1,1000);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
