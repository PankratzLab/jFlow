package org.genvisis.one.george;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashMap;

public class HumanLabels {
	// filename is .csv containing human labels for each data point

	// want cumulative, not just leaf
	// val: lymphocytes, single lymphs, live single lymphs, B cells, T cells, cytotoxic T cells,  IgD+ mem Bcells  
	
	ArrayList<Integer[]> data = new ArrayList<Integer[]>();
	public HumanLabels(String filepath) {
		Node[] helperTChildren = new Node[]{new Node(17, null), new Node(27, 38, null), new Node(20, 25, null), new Node(19, 37, null), new Node(34, null)};
		Node helperT = new Node(10, helperTChildren);
		
		Node[] cytotoxicTEffectorChildren = new Node[]{new Node(29, null), new Node(24, null), new Node(33, null)};
		Node[] cytotoxicTEffectorMemoryChildren = new Node[]{new Node(32, null), new Node(22, null), new Node(28, null), new Node(23, null)};
		Node[] cytotoxicTChildren = new Node[]{new Node(26, 30, null), new Node(16, cytotoxicTEffectorChildren), new Node(15, 31, null), new Node(18, 35, cytotoxicTEffectorMemoryChildren)};

		Node cytotoxicT = new Node(7, cytotoxicTChildren);
		Node tCellNode = new Node(6, new Node[]{cytotoxicT, helperT});
		
		Node[] bCellChildren= new Node[]{new Node(11, null), new Node(8, null), new Node(9, null)};
		Node bCellNode = new Node(5, bCellChildren);
		
		Node liveSingleNode = new Node(3, new Node[]{tCellNode, bCellNode});
		Node singleNode = new Node(2, new Node[]{liveSingleNode});
		Node lymphocyteNode = new Node(1, new Node[]{singleNode});

		Tree panel1Tree = new Tree(lymphocyteNode, 38);

		File f = new File(filepath);
		try {
			Scanner s = new Scanner(f);
			s.useDelimiter("\n");
			
			while (s.hasNextInt()) {
				Integer d = s.nextInt();
				try {
					Integer[] arr = panel1Tree.getArray(d);
					// LEAFID NOT FOUND
					if (arr == null) {
						System.out.println("LEAFID NOT FOUND setting to 0's: " + d);
						arr = new Integer[panel1Tree.maxLeafId];
						for (int i = 0; i < arr.length; i++)
							arr[i] = 0;
					}
					data.add(arr);
				}
				catch (Exception e) {
					System.out.println("err: leaf id not in tree");
				}
			}
			s.close();
		}
		catch (FileNotFoundException e) {
			System.out.println("err: file not found");
		}
	}
	
	public void exportData(String out_filename) {
		try {
			PrintWriter w = new PrintWriter(out_filename, "UTF-8");
			StringBuilder out = new StringBuilder();
			for (int i = 0; i < data.size(); i++) {
				Integer[] arr = data.get(i);
				for (int j = 0; j < arr.length; j++) {
					int v = arr[j];
					out.append(new Integer(v).toString());

					if (j != arr.length-1) {
						out.append(", ");
					}
					else {
						out.append("\n");
					}
				}				
			}
			w.write(out.toString());
			w.close();
		}
		catch (IOException e) {
			System.out.println("err: invalid filename");
		}
	}
}
