package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class FixRotterdam {
	public static void main(String[] args) {
		BufferedReader reader;
		PrintWriter writer;

		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Folson\\METAL\\Rotterdam_src\\";
		String file1 = "RS2VTE_CHARGE.chargefmt.RS.v3.txt";
		String file2 = "RS2VTE_CHARGE.chargefmt.RS.v2.txt";
		String outfile = "RS2VTE_CHARGE.chargefmt.RS.v4.txt";
		try {
			reader = new BufferedReader(new FileReader(dir + file1));
			writer = new PrintWriter(new FileWriter(dir + outfile));
			for (int i = 0; i < 414388; i++) {
				writer.println(reader.readLine());
			}
			reader.close();

			reader = new BufferedReader(new FileReader(dir + file2));
			for (int i = 0; i < 193555; i++) {
				reader.readLine();
			}
			while (reader.ready()) {
				writer.println(reader.readLine());
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + file1 + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + file1 + "\"");
			System.exit(2);
		}

	}
}
