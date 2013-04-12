package filesys;

import java.io.*;

import common.*;

public class FamilyStructure {
	public static final String[][] TYPICAL_HEADERS = {{"FID", "famid"}, {"IID", "id"}, {"fa"}, {"mo"}, {"sex"}}; 
	
	private String[][] ids;
	private byte[] genders;
	private byte[] affections;
	private String[] dnas;

	public FamilyStructure(String[][] ids, byte[] genders, byte[] affections) {
		this(ids, genders, affections, null);
	}
	
	public FamilyStructure(String[][] ids, byte[] genders, byte[] affections, String[] dnas) {
		this.ids = ids;
		this.genders = genders;
		this.affections = affections;
		this.dnas = dnas;
	}
	
	public FamilyStructure(String filename) {
		this(filename, false);
	}
	
	public FamilyStructure(String filename, boolean loadDNAs) {
		BufferedReader reader;
		String[] line;
		int count;

		try {
			reader = new BufferedReader(new FileReader(filename));
			count = 0;
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();

			reader = new BufferedReader(new FileReader(filename));
			ids = new String[count][];
			genders = new byte[count];
			affections = new byte[count];
			dnas = loadDNAs?new String[count]:null;
			for (int i = 0; i<count; i++) {
				line = reader.readLine().trim().split("[\\s]+");
				ids[i] = new String[] {line[0], line[1], line[2], line[3]};
				genders[i] = Byte.parseByte(line[4]);
				affections[i] = Byte.parseByte(line[5]);
				if (loadDNAs) {
					dnas[i] = line[6];
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}
	}
	
	public String[][] getIds() {
		return ids;
	}

	public byte[] getGenders() {
		return genders;
	}

	public byte[] getAffections() {
		return affections;
	}

	public void setAffections(byte[] affs) {
		affections = affs;
	}
	
	public String[] getDnas() {
		return dnas;
	}
	
	public String getIndividualHeader(int index, boolean displayDNA) {
		return Array.toStr(ids[index])+"\t"+genders[index]+"\t"+affections[index]+(displayDNA&&dnas == null?"\t"+dnas[index]:"");
	}

	public void writeToFile(String filename, boolean displayDNA) {
        PrintWriter writer;
        
        try {
	        writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i<ids.length; i++) {
				writer.println(Array.toStr(ids[i])+"\t"+genders[i]+"\t"+affections[i]+(displayDNA&&dnas == null?"\t"+dnas[i]:""));
	        }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }
	}
	
	public static boolean likelyPedHeader(String[] line) {
		return Array.countIf(ext.indexFactors(TYPICAL_HEADERS, line, false, true, false, false), -1) < 3;
	}
}
