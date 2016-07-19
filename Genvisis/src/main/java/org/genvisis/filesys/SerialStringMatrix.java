// -Xms1024M -Xmx1024M

package filesys;

import java.io.Serializable;
import java.io.*;
import common.*;

public class SerialStringMatrix implements Serializable {
	public static final long serialVersionUID = 1L;
	private String[][] matrix;

	public String[][] getMatrix() {
		return matrix;
	}
	
	public SerialStringMatrix(String[][] matrix) {
		this.matrix = matrix;
	}
	
	public SerialStringMatrix(String filename, String delimiter) {
		BufferedReader reader;
        int count;
        
        try {
        	count = Files.countLines(filename, 0);
        	matrix = new String[count][];
	        reader = new BufferedReader(new FileReader(filename));
	        count = 0;
	        while (reader.ready()) {
	        	matrix[count] = reader.readLine().trim().split(delimiter);
	        	count++;
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

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SerialStringMatrix load(String filename, boolean jar) {
		return (SerialStringMatrix)Files.readSerial(filename, jar, true);
	}
 	
	public static void main(String[] args) {
		String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\ARIC_autosomes_table.csv";
		new SerialStringMatrix(filename, ",").serialize(filename+".ser");
    }
}
