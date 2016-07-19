package link.bat;

import java.io.*;
import java.util.*;

public class dat2loc {

	public dat2loc() {}

	public dat2loc(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String temp;
		StringTokenizer st;
		int numMarkers, numAlleles;
		float[] mrkrDists;
		float fTemp, multiplier;

		try {
			reader = new BufferedReader(new FileReader(filename));
		} catch (Exception e) {
			try {
				reader = new BufferedReader(new FileReader("/home/npankrat/park/00masters/"+filename));
				System.err.println("Could not find "+filename+" in the current directory");
				System.err.println("  using the one in /home/npankrat/park/00masters/");
			} catch (Exception e2) {
				System.err.println("Could not find "+filename+" in /home/npankrat/park/00masters/ or in the current directory");
				System.exit(1);
			}

		}

		if (filename.startsWith("map")) {
			writer = new PrintWriter(new FileWriter(filename.substring(0, filename.length()-4)+".loc"));
		} else if (filename.startsWith("nuke")) {
			writer = new PrintWriter(new FileWriter("map"+filename.substring(filename.length()-2)+".loc"));
		} else {
			writer = new PrintWriter(new FileWriter(filename+".loc"));
		}

		st = new StringTokenizer(reader.readLine());

		numMarkers = Integer.valueOf(st.nextToken()).intValue()-1;
		writer.println(numMarkers);
		writer.println();

		reader.readLine();
		reader.readLine();
		temp = reader.readLine();
		st = new StringTokenizer(temp);
		temp = st.nextToken();
		if (temp.equals("1")) { // if it is a qualitative phenotype, then pass 3
			// more
			for (int i = 0; i<3; i++) {
				reader.readLine();
			}
			if (filename.indexOf("23")!=-1) {
				reader.readLine();
			}
		}
		if (temp.equals("0")) { // if it is a quantative phenotype, then pass 5
			// more
			for (int i = 0; i<5; i++) {
				reader.readLine();
			}
		}
		for (int i = 0; i<numMarkers; i++) {
			reader.readLine();
			reader.readLine();
		}
		reader.readLine();

		st = new StringTokenizer(reader.readLine());
		fTemp = Float.valueOf(st.nextToken()).floatValue();
		if (fTemp<1) {
			multiplier = 100;
		} else {
			multiplier = 1;
		}
		mrkrDists = new float[numMarkers];
		for (int i = 0; i<numMarkers-1; i++) {
			mrkrDists[i] = Float.valueOf(st.nextToken()).floatValue()*multiplier;
			if (mrkrDists[i]<0.1) {
				mrkrDists[i] = (float)0.1;
			}
		}

		reader.close();
		reader = new BufferedReader(new FileReader(filename));

		reader.readLine();
		reader.readLine();
		reader.readLine();
		temp = reader.readLine();
		st = new StringTokenizer(temp);
		temp = st.nextToken();
		if (temp.equals("1")) { // if it is a qualitative phenotype, then pass 3
			// more
			for (int i = 0; i<3; i++) {
				reader.readLine();
			}
			if (filename.indexOf("23")!=-1) {
				reader.readLine();
			}
		}
		if (temp.equals("0")) { // if it is a quantative phenotype, then pass 5
			// more
			for (int i = 0; i<5; i++) {
				reader.readLine();
			}
		}

		for (int i = 0; i<numMarkers; i++) {
			st = new StringTokenizer(reader.readLine());
			st.nextToken();
			numAlleles = Integer.valueOf(st.nextToken()).intValue();
			st.nextToken();
			writer.println(st.nextToken()+" "+numAlleles);
			for (int j = 1; j<=numAlleles; j++) {
				writer.print(j+" ");
			}
			writer.println();
			st = new StringTokenizer(reader.readLine(), "<");
			writer.println(st.nextToken());
			if (i!=numMarkers-1) {
				writer.println();
				writer.println((mrkrDists[i]));
				writer.println();
			}
		}
		reader.close();
		writer.close();

	}

	public static void main(String[] args) throws IOException {
		try {
			new dat2loc(args[0]);
			// stdin.readLine();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
