package stats;

import java.util.ArrayList;
//import java.util.Hashtable;
//import java.util.Vector;

import common.AlleleFreq;
import common.Array;
import common.ext;
//import common.CountVector;
//import common.HashVec;
//import common.Sort;

public class CTable {
	/*
	private String[] rowLabels;
	private String[] columnLabels;
	private int[][] cells;
	*/
	private ArrayList<String> rowLabels;
	private ArrayList<String> columnLabels;
	private ArrayList<ArrayList<Integer>> cells;
	private int[][] cellsCondensedForCallrate;
	int numberOfNAColumns;
	
	public CTable() {}
	
//	public CTable(int[][] matrix) {
//		this(matrix[0], matrix[1]);
//	}
	

	public CTable(int[] array1, int[] array2) {
		this(Array.toStringArray(array1), Array.toStringArray(array2));
	}

	public CTable(int[] array1, int[] array2, String[][] array1ClassLabels, String[][] array2ClassLabels) {
		this(array1, array2);
		replaceIdWithLabel(array1ClassLabels, array2ClassLabels);
	}

	public CTable(String[] array1, String[] array2) {
		rowLabels = new ArrayList<String>();
		columnLabels = new ArrayList<String>();
		cells = new ArrayList <ArrayList<Integer>>();
		
		int rowIndex;
		int columnIndex;
		
		if (array1.length != array2.length) {
			System.err.println("Error - arrays are not of the same length");
		}
		
		for(int i=0; i<array1.length; i++) {
			rowIndex=rowLabels.indexOf(array1[i]);
			if (rowIndex<0) {
				rowIndex=0;
				for (int j=0; j<rowLabels.size(); j++){
					if(Integer.parseInt(array1[i])<Integer.parseInt(rowLabels.get(j))){
						rowIndex=j;
						break;
					} else {
						rowIndex=j+1;
					}
				}
				rowLabels.add(rowIndex,array1[i]); //to put sorting here. But need to adjust cells[][] as well
				//rowIndex=rowLabels.size()-1;
				cells.add(rowIndex, new ArrayList<Integer>());
				for (int j=0; j<columnLabels.size(); j++) {
					cells.get(rowIndex).add(0);
				}
			}
			columnIndex=columnLabels.indexOf(array2[i]);
			if (columnIndex<0) {
				columnIndex=0;
				for (int j=0; j<columnLabels.size(); j++) {
					if (Integer.parseInt(array2[i])<Integer.parseInt(columnLabels.get(j))) {
						columnIndex=j;
						break;
					} else {
						columnIndex=j+1;
					}
				}
				columnLabels.add(columnIndex, array2[i]); //to put sorting here. But need to adjust cells[][] as well
				for (int j=0; j<cells.size(); j++) {
					cells.get(j).add(columnIndex,0);
				}
			}
			cells.get(rowIndex).set(columnIndex, cells.get(rowIndex).get(columnIndex)+1);
			//System.out.println("array1: "+array1[i]+"\t array2: "+array2[i]+"\t rowLabels size: "+rowLabels.size()+"\t columnLabels size: "+columnLabels.size()+"\t cells size: "+cells.size()+" x "+cells.get(0).size()+"\t rowIndex: "+rowIndex+"\t columnIndex: "+columnIndex);
			/*
			System.out.println("array1: "+array1[i]+"\t array2: "+array2[i]+"\t cells["+rowIndex+"]["+columnIndex+"]: "+cells.get(rowIndex).get(columnIndex));
			for (int x=0; x<cells.size(); x++) {
				for (int y=0; y<cells.get(0).size(); y++) {
					System.out.println("cells["+x+"]["+y+"]: "+cells.get(x).get(y));
				}
			}
			*/
		}
	}
	
	/*
	public int[][] createTable(boolean treatAsNumbers) {
		String[][] rowLookup;
		String[][] columnLookup;
		String[] sortedRows;
		String[] sortedColumns;

		
		
		sortedRows = Sort.putInOrder(Array.toStringArray(rowLabels), treatAsNumbers);
		
		for (int i=0; i<sortedRows.length; i++) {
			rowLookup[i][0] = rowLabels.indexOf(sortedRows[i]);
			rowLookup[i][1] = rowLabels.indexOf(sortedRows[i]);
		}
		
        createTable
	}
	*/
	
	public void replaceIdWithLabel(String[][] rowLookup, String[][] columnLookup) {
		for (int i=0; i<rowLabels.size(); i++) {
			for (int j=0; j<rowLookup.length; j++) {
				//System.out.println("Looking up rowLabel. rowLabels: "+rowLabels.get(i)+"\t rowLookup: "+rowLookup[j][0]+"\t Equal? "+(rowLabels.get(i).contentEquals(rowLookup[j][0])));
				if (rowLabels.get(i).contentEquals(rowLookup[j][0])) {
					rowLabels.set(i, rowLookup[j][1]);
					break;
				}
				if (j==rowLookup.length-1) {
					rowLabels.set(i, "N/A"); //if the rowID does not have a rowLabel value, the set the rowID's lable to N/A.
				}
			}
		}
		for (int i=0; i<columnLabels.size(); i++) {
			for (int j=0; j<columnLookup.length; j++) {
				if (columnLabels.get(i).contentEquals(columnLookup[j][0])) {
					columnLabels.set(i, columnLookup[j][1]);
					break;
				}
				if (j==columnLookup.length-1) {
					columnLabels.set(i, "N/A");
				}
			}
		}
	}
	
	public String generateToolTipText() {
		String output;
		output = "<html><table border=\"1\">";// <tr><td></td>" + columnLabels.get(1) + "</tr>";
		for (int i=-1; i<rowLabels.size(); i++) {
			for (int j=-1; j<columnLabels.size(); j++) {
				if (i==-1) {
					if (j==-1) {
						output = output + "<tr><td></td>";
					} else {
						output = output + "<td>" + columnLabels.get(j) + "</td>";
					}
				} else if (j==-1) {
					output = output + "<tr><td align=\"center\">"+rowLabels.get(i) + "</td>";
				} else {
					output = output + "<td align=\"center\">" + cells.get(i).get(j) + "</td>";
				}
			}
			output = output + "</tr>";
		}
		output = output + "</table>" + "</html>";
		return output;
	}

	public String generateToolTipTextForCallRate() {
		String output;
		output = "<html><table border=\"1\">";// <tr><td></td>" + columnLabels.get(1) + "</tr>";
		for (int i=-1; i<2; i++) {
			for (int j=-1; j<cellsCondensedForCallrate[0].length; j++) {
				if (i==-1) {
					if (j==-1) {
						output = output + "<tr><td></td>";
					} else {
						output = output + "<td>" + columnLabels.get(j + numberOfNAColumns) + "</td>";
					}
				} else if (j==-1) {
					output = output + "<tr><td align=\"center\">" + (i==0?"Missing":"A/A, A/B, or B/B") + "</td>";
				} else {
					output = output + "<td align=\"center\">" + cellsCondensedForCallrate[i][j] + "</td>";
				}
			}
			output = output + "</tr>";
		}
		output = output + "</table>" + "</html>";
		return output;
	}

	public int[][] getContingencyTable() {
		int[][] output = new int[rowLabels.size()][columnLabels.size()];
		for (int i=0; i<rowLabels.size(); i++) {
			for (int j=0; j<columnLabels.size(); j++) {
				output[i][j]= cells.get(i).get(j);
			}
		}
		return output;
	}
	
	public int[][] getContingencyTableForCallRate() {
		numberOfNAColumns = 0;
		for (int i=0; i<columnLabels.size(); i++) {
			if (columnLabels.get(i).contentEquals("N/A")) {
				numberOfNAColumns ++;
			}
		}
		cellsCondensedForCallrate = new int[2][columnLabels.size()-numberOfNAColumns];
		for (int i=0; i<rowLabels.size(); i++) {
			for (int j=0; j<columnLabels.size(); j++) {
				if (columnLabels.get(j).contentEquals("N/A")) {
				} else if (rowLabels.get(i).contentEquals("N/A")) {
					cellsCondensedForCallrate[0][j-numberOfNAColumns]= cellsCondensedForCallrate[0][j-numberOfNAColumns] + cells.get(i).get(j);
				} else {
					cellsCondensedForCallrate[1][j-numberOfNAColumns]= cellsCondensedForCallrate[1][j-numberOfNAColumns] + cells.get(i).get(j);
				}
			}
		}
		return cellsCondensedForCallrate;
	}
	
	public float getMinorAlleleFrequency() {
		float result = 0;
		int[] AB = new int[3];
		int x =0;
		for (int i=0; i<rowLabels.size(); i++) {
			if (rowLabels.get(i).contentEquals("N/A")) {
				x=1;
			} else {
				for (int j=0; j<columnLabels.size(); j++) {
					if (columnLabels.get(j).contentEquals("N/A")) {
					} else {
						AB[i-x]= AB[i-x] + cells.get(i).get(j);
					}
				}
			}
		}
		result = (float)(2*AB[0]+AB[1])/(float)(2*(AB[0]+AB[1]+AB[2]));
		if (result>0.5) {
			result = 1 - result;
		}
		return result;
	}

	
	public int[][] getAlleleFreqBySex() {
		int[][] result = new int[2][2];
		int rowIndex, columnIndex;
		if (rowLabels.get(0).contentEquals("N/A")) {
			rowIndex=1;
		} else {
			rowIndex=0;
		}
		if (columnLabels.get(0).contentEquals("N/A")) {
			columnIndex=1;
		} else {
			columnIndex=0;
		}
		//result = new int[][] {{cells.get(rowIndex).get(columnIndex+1), cells.get(rowIndex+1).get(columnIndex+1), cells.get(rowIndex+2).get(columnIndex+1)},
		//						{cells.get(rowIndex).get(columnIndex), cells.get(rowIndex+1).get(columnIndex), cells.get(rowIndex+2).get(columnIndex)}};
		result[0][0] = cells.get(rowIndex).get(columnIndex+1)*2 + cells.get(rowIndex+1).get(columnIndex+1);
		result[0][1] = cells.get(rowIndex+1).get(columnIndex+1) + cells.get(rowIndex+2).get(columnIndex+1)*2;
		result[1][0] = cells.get(rowIndex).get(columnIndex)*2 + cells.get(rowIndex+1).get(columnIndex);
		result[1][1] = cells.get(rowIndex+1).get(columnIndex) + cells.get(rowIndex+2).get(columnIndex)*2;
		return result;
	}
	
	public String generateToolTipTextForAllelFreq() {
		String output;
		int rowIndex, columnIndex;
		if (rowLabels.get(0).contentEquals("N/A")) {
			rowIndex=1;
		} else {
			rowIndex=0;
		}
		if (columnLabels.get(0).contentEquals("N/A")) {
			columnIndex=1;
		} else {
			columnIndex=0;
		}
		output = "<html><table border=\"1\">" +
					"<tr>" +
						"<td></td>" +
						"<td align=\"center\"> Female  </td>" +
						"<td align=\"center\"> Male </td>" +
					"</tr>" +
					"<tr>" +
						"<td align=\"center\"> Allele Frequency </td>" +
						"<td align=\"center\">"+ext.formDeci(AlleleFreq.calcFrequency(cells.get(rowIndex).get(columnIndex+1), cells.get(rowIndex+1).get(columnIndex+1), cells.get(rowIndex+2).get(columnIndex+1)),4)+"</td>" +
						"<td align=\"center\">"+ext.formDeci(AlleleFreq.calcFrequency(cells.get(rowIndex).get(columnIndex), cells.get(rowIndex+1).get(columnIndex), cells.get(rowIndex+2).get(columnIndex)),4)+"</td>" +
					"</tr>";
		return output;
	}



	/*
	public CTable(int[] array1, int[] array2) {
		//Hashtable<String,String> hash;
		ArrayList v1, v2 = new ArrayList();
		
		cells = new int[1][1];

		if (array1.length != array2.length) {
			System.err.println("Error - arrays are not of the same length");
		}
		
		for(int i=0; i<array1.length; i++) {
			if (array1[i]<=v1.size() || array2[i]<=v2.size()) {
				v1.add(array1[i]);
				v2.add(array2[i]);
			}
			cells[array1[i]][array2[i]] ++;
		}
		rowLabels = Sort.putInOrder(Array.toStringArray(v1));
		columnLabels = Sort.putInOrder(Array.toStringArray(v2));
	}
	*/

	
	/*
	public CTable(String[] array1, String[] array2) {
		Hashtable<String,String> hash;
		Vector<String> v1, v2;
		String trav;
		int count;
		
		hash = new Hashtable<String, String>();

		if (array1.length != array2.length) {
			System.err.println("Error - arrays are not of the same length");
		}
		
		v1 = new Vector<String>();
		v2 = new Vector<String>();
		for(int i=0; i<array1.length; i++) {
			trav = array1[i]+" "+array2[i];
			if (hash.containsKey(trav)) {
				count = Integer.parseInt(hash.get(trav));
			} else {
				count = 0;
			}
			count++;
			hash.put(trav, count+"");
			
			HashVec.addIfAbsent(array1[i], v1);
			HashVec.addIfAbsent(array2[i], v2);
		}
		rowLabels = Sort.putInOrder(Array.toStringArray(v1));
		columnLabels = Sort.putInOrder(Array.toStringArray(v2));
		cells = new int[1][1];//???
	}
	*/


}
