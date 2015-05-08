package one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import cnv.manage.PlinkData;
import common.Elision;
import common.Files;
import common.HashVec;
import common.ext;

public class PlinkMarkerLoader {
	
	public static void main(String[] args) {
		String[] markerNames = HashVec.loadFileToStringArray("D:/PlinkGeno/mkrs10000.txt", false, null, false);
		String plinkFileRoot = "D:/PlinkGeno/plink";
		byte[][] genotypes = (new PlinkMarkerLoader()).run(plinkFileRoot, markerNames);
		System.out.println(ext.getTime() + "]\tFinished");
	}
	
	public static int[] markerLookup(String plinkFileRoot, String[] markers) {
		BufferedReader reader;
		String[] line;
		HashSet<String> lookFor = new HashSet<String>();
		for (String s : markers) {
			lookFor.add(s);
		}
		HashMap<String, Integer> markerIndicesTemp = new HashMap<String, Integer>();
		int[] markerIndices = new int[markers.length];
		int cnt = 0;
		
		try {
			reader = new BufferedReader(new FileReader(plinkFileRoot+".bim"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\s+");
				String mkr = line[1];
				if (lookFor.contains(mkr)) {
					markerIndicesTemp.put(mkr, cnt);
				}
				cnt++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + plinkFileRoot+".bim" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + plinkFileRoot+".bim" + "\"");
			System.exit(2);
		}
		
		for (int i = 0; i < markers.length; i++) {
			markerIndices[i] = markerIndicesTemp.get(markers[i]).intValue();
		}
		
		return markerIndices;
	}

	private String[] sortMarkers(String[] markers, int[] positions) {
		final HashMap<String, Integer> unsortedMap = new HashMap<String, Integer>();
		for (int i = 0; i < markers.length; i++) {
			unsortedMap.put(markers[i], positions[i]);
		}
		TreeMap<String, Integer> sortedMap = new TreeMap<String, Integer>(new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				return unsortedMap.get(o1).compareTo(unsortedMap.get(o2));
			}
		});
		sortedMap.putAll(unsortedMap);
		String[] sortedMarkers = new String[markers.length];
		Set<String> sortedKeys = sortedMap.keySet();
		return sortedKeys.toArray(sortedMarkers);
	}
	
	private int[] sortMarkerIndices(int[] positions) {
		int[] sortedPositions = new int[positions.length];
		ArrayList<Integer> sort = new ArrayList<Integer>();
		for (int pos : positions) {
			sort.add(pos);
		}
		sort.sort(new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return o1.compareTo(o2);
			}
		});
		for (int i = 0; i < sort.size(); i++) {
			sortedPositions[i] = sort.get(i).intValue();
		}
		return sortedPositions;
	}
	
	private byte[][] run(String plinkDirAndFilenameRoot, String[] markers) {
		RandomAccessFile in = null;

		HashMap<String, byte[]> mkrGenotypes = new HashMap<String, byte[]>();
		byte[][] genoBytes = new byte[markers.length][];
		
		try {
			in = new RandomAccessFile(plinkDirAndFilenameRoot + ".bed", "r");

			byte[] magicBytes = new byte[3];
			in.read(magicBytes);
			if (magicBytes[2] == 0) {
				// sample dominant
				System.err.println("Error - .bed file is sample-dominant.");
			} else {
				System.out.println(ext.getTime() + "]\tLoading sample count");
				int famCnt = Files.countLines(plinkDirAndFilenameRoot + ".fam", 0); 
				int blockSize = (int) ((double) famCnt / 4.0d); // TODO check for non-completeness (i.e. N not evenly divisible by 4)
				System.out.println(ext.getTime() + "]\tLoading marker count");
//				int V = Files.countLines(plinkDirAndFilenameRoot + ".bim", 0); // # of markers
				System.out.println(ext.getTime() + "]\tLoading positions for given markers");
				int[] markerIndicesUnsorted = markerLookup(plinkDirAndFilenameRoot, markers);
				System.out.println(ext.getTime() + "]\tSorting given markers by position");
				String[] markersSorted = sortMarkers(markers, markerIndicesUnsorted);
				int[] markerIndices = sortMarkerIndices(markerIndicesUnsorted);
				
//				int totalBytes = V * blockSize;
				
				for (int i = 0; i < markerIndices.length; i++) {
					in.seek(markerIndices[i] * blockSize);
					
					byte[] markerBytes = new byte[blockSize];
					byte[] sampGeno = new byte[famCnt];
					in.read(markerBytes);
					for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
						byte bedByte = markerBytes[bitInd];
						byte[] genotypes = PlinkData.decodeBedByte(bedByte);
						for (int g = 0; g < genotypes.length; g++) {
							sampGeno[bitInd * 4 + g] = genotypes[g];
						}
					}
					
					mkrGenotypes.put(markersSorted[i], sampGeno);
					
				}
			}
			
			in.close();
			
			for (int i = 0; i < markers.length; i++) {
				genoBytes[i] = mkrGenotypes.get(markers[i]);
			}
			
			return genoBytes;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Elision e) {
			e.printStackTrace();
		} finally {
			if (in != null) {
				try {
					in.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		return null;
		
	}
	
	
}
