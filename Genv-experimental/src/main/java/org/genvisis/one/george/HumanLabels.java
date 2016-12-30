package org.genvisis.one.george;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.HashMap;

public class HumanLabels {
	// filename is .csv containing human labels for each data point
	HashMap<Integer, Integer[]> mapping = new HashMap<Integer, Integer[]>();
	// want cumulative, not just leaf
	// val: lymphocytes, single lymphs, live single lymphs, B cells, T cells, cytotoxic T cells,  IgD+ mem Bcells  
	
	public HumanLabels(String filename) {
		Node[] helperTChildren = new Node[]{new Node(17, null), new Node(27, 38, null), new Node(20, 25, null), new Node(19, 37, null), new Node(34, null)};
		Node helperT = new Node(10, helperTChildren);
		
		//Node[] cytotoxicTChildren=
		//Node cytotoxicT
		
		//Node[] bCellChildren=
		//Node bCellNode
		
		//Node tCellNode = 
		//Node liveSingleNode =
		//Node singleNode =
		//Node lymphNode =
		
		File f = new File(filename);
		try {
			Scanner s = new Scanner(f);
			s.useDelimiter("\n");
			
			String outStr = "";
			while (s.hasNextInt()) {
				Integer d = s.nextInt();
				outStr = outStr + d.toString() + " ";
			}
			s.close();
			System.out.print(outStr);
		}
		catch (FileNotFoundException e) {
			System.out.println("err: file not found");
		}
	}
	
	public void exportToMatrix(String filename) {
		FileWriter w = FileWriter(filename);
		w.write()
	}
}
