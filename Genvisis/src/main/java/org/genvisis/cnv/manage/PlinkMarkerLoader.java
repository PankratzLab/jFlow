package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class PlinkMarkerLoader implements Runnable {

	public static void main(String[] args) {
		String[] markerNames = HashVec.loadFileToStringArray(	"D:/PlinkGeno/mkrs10000.txt", false, null, false);
		String plinkFileRoot = "D:/PlinkGeno/plink";
		(new PlinkMarkerLoader(null, plinkFileRoot, markerNames)).run();
		System.out.println(ext.getTime() + "]\tFinished");
	}

	String fileRoot;
	String[] markerList;
	String[] famIDList;
	int[] markerPosInBim;
	volatile byte[][] genotypes;
	volatile boolean[] loaded;
	Logger log;
	HashMap<String, Integer> markerIndicesLookup;
	HashMap<String, Integer> famIDLookup;

	boolean idListsDiffer;
	boolean initialized;
	boolean killed;
	boolean killComplete;
	Thread thread;

	public PlinkMarkerLoader(Project proj, String plinkFileRoot, String[] markers) {
		fileRoot = plinkFileRoot;
		markerList = markers;
		log = proj == null ? new Logger() : proj.getLog();

		if (markers == null) {
			log.reportError("The list of markers for MarkerDataLoader to load was null");
			killed = true;
			return;
		}

		if (markers.length == 0) {
			log.reportError("The list of markers for MarkerDataLoader to load was empty (n=0)");
			killed = true;
			return;
		}

		markerPosInBim = new int[markerList.length];
		genotypes = new byte[markerList.length][];
		loaded = ArrayUtils.booleanArray(markerList.length, false);

		lookupMarkerPositions();
		lookupIDs();
	}

	private void initiate() {
		initialized = true;
	}

	public boolean isKilled() {
		return killed;
	}

	private void registerThread(Thread thread) {
		this.thread = thread;
	}

	public void kill() {
		killed = true;
		if (!thread.isAlive()) {
			killComplete = true;
		}
	}

	public boolean killComplete() {
		return killComplete;
	}

	public Thread getThread() {
		return thread;
	}

	private void lookupMarkerPositions() {
		BufferedReader reader;
		String[] line;
		HashSet<String> lookFor = new HashSet<String>();
		for (String s : markerList) {
			lookFor.add(s);
		}
		markerIndicesLookup = new HashMap<String, Integer>();
		HashMap<String, Integer> markerIndicesLookupTemp = new HashMap<String, Integer>();
		int cnt = 0;
		try {
			reader = Files.getAppropriateReader(fileRoot + ".bim");
			String temp;
			while ((temp = reader.readLine()) != null) {
				line = temp.trim().split("\\s+");
				String mkr = line[1];
				if (lookFor.contains(mkr)) {
					markerIndicesLookupTemp.put(mkr, cnt);
				}
				cnt++;
			}
			reader.close();

			for (int i = 0; i < markerList.length; i++) {
				markerPosInBim[i] = markerIndicesLookupTemp.get(markerList[i]) == null ? -1 : markerIndicesLookupTemp.get(markerList[i]).intValue();
				markerIndicesLookup.put(markerList[i], i);
			}
		} catch (FileNotFoundException fnfe) {
			// TODO should KILL here
			log.reportException(fnfe);
		} catch (IOException ioe) {
			// TODO should KILL here
			log.reportException(ioe);
		}
	}

	private void lookupIDs() {
		BufferedReader reader;
		String[] line;

		famIDList = new String[Files.countLines(fileRoot + ".fam", 0)];

		try {
			reader = Files.getAppropriateReader(fileRoot + ".fam");
			String temp = null;
			int cnt = 0;
			while ((temp = reader.readLine()) != null) {
				line = temp.trim().split("[\\s]+");
				famIDList[cnt] = line[0] + "\t" + line[1];
				cnt++;
			}
			reader.close();
			reader = null;
		} catch (FileNotFoundException fnfe) {
			// TODO should KILL here
			log.reportException(fnfe);
		} catch (IOException ioe) {
			// TODO should KILL here
			log.reportException(ioe);
		}

		famIDLookup = new HashMap<String, Integer>();
		for (int i = 0; i < famIDList.length; i++) {
			famIDLookup.put(famIDList[i], i);
		}
	}

	@Override
	public void run() {
		RandomAccessFile in = null;

		HashMap<String, byte[]> mkrGenotypes = new HashMap<String, byte[]>();

		if (killed) {
			return;
		}

		initiate();
		int cnt = 0;
		long time = new Date().getTime();

		try {
			in = new RandomAccessFile(fileRoot + ".bed", "r");

			byte[] magicBytes = new byte[3];
			in.read(magicBytes);
			if (magicBytes[2] == 0) {
				log.reportError("Error - .bed file is sample-dominant.");
			} else {
				int famCnt = Files.countLines(fileRoot + ".fam", 0);
				int blockSize = (int) Math.ceil(famCnt / 4.0d);
				for (int i = 0; i < markerList.length; i++) {
					if (markerPosInBim[i] == -1) {
						// missing marker, not present in PLINK files
						mkrGenotypes.put(markerList[i], ArrayUtils.byteArray(famIDList.length, (byte) -1));
						cnt++;
						continue;
					}
					in.seek(3 + markerPosInBim[i] * blockSize);

					byte[] markerBytes = new byte[blockSize];
					byte[] sampGeno = new byte[famIDList.length];
					in.read(markerBytes);
					for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
						byte bedByte = markerBytes[bitInd];
						byte[] genos = PlinkData.decodeBedByte(bedByte);

						for (int g = 0; g < genos.length; g++) {
							int idInd = (bitInd * 4 + g) >= famIDList.length ? -1 : famIDLookup.get(famIDList[bitInd * 4 + g]) == null ? -1 : famIDLookup.get(famIDList[bitInd * 4 + g]);
							if (idInd == -1 || idInd > sampGeno.length) {
								continue;
							}
							sampGeno[idInd] = genos[g];
						}
					}

					mkrGenotypes.put(markerList[i], sampGeno);
					cnt++;
				}
			}

			in.close();

			for (int i = 0; i < markerList.length; i++) {
				genotypes[i] = mkrGenotypes.get(markerList[i]);
				loaded[i] = true;
			}

		} catch (FileNotFoundException e) {
			log.reportException(e);
			// TODO should KILL here
		} catch (IOException e) {
			log.reportException(e);
			// TODO should KILL here
		} catch (Elision e) {
			log.reportException(e);
			// TODO should KILL here
		} finally {
			if (in != null) {
				try {
					in.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		int aa, ab, bb, m;
		aa = 0;
		ab = 0;
		bb = 0;
		m = 0;
		for (int i = 0; i < genotypes[0].length; i++) {
			byte g = genotypes[0][i];
			if (g == 0) {
				aa++;
			} else if (g == 1) {
				ab++;
			} else if (g == 2) {
				bb++;
			} else if (g == -1) {
				m++;
			}
		}
		System.out.println("PLINK GENOTYPES  --->  AA: " + aa + " AB: " + ab + " BB: " + bb + " Miss: " + m);

		if (killed) {
			log.report("PlinkMarkerLoader killed");
			markerList = null;
			markerPosInBim = null;
			loaded = null;
			genotypes = null;
			System.gc();
			killComplete = true;
		} else {
			log.report("Independent thread has finished loading "	+ cnt + " markers in "
									+ ext.getTimeElapsed(time));
		}

	}

	public byte getGenotypeForIndi(String marker, String fidiid) {
		int idIndex = famIDLookup.get(fidiid) == null ? -1 : famIDLookup.get(fidiid);
		if (idIndex == -1) {
			return (byte) -1;
		} else {
			int markerIndex = markerIndicesLookup.get(marker) == null	? -1 : markerIndicesLookup.get(marker);
			if (markerIndex == -1) {
				return (byte) -1;
			}
			while (!loaded[markerIndex]) {
				Thread.yield();
			}
			return genotypes[markerIndex][idIndex]; 
		}
	}

	public static PlinkMarkerLoader loadPlinkDataFromListInSeparateThread(Project proj,
																																				String plinkDirFileRoot,
																																				String[] markerList) {
		PlinkMarkerLoader plinkMarkerLoader;
		Thread thread;

		proj.getLog().report("PLINK marker data is loading in an independent thread.");
		plinkMarkerLoader = new PlinkMarkerLoader(proj, plinkDirFileRoot, markerList);
		if (plinkMarkerLoader.isKilled()) {
			return null;
		}
		plinkMarkerLoader.initiate();
		thread = new Thread(plinkMarkerLoader);
		thread.start();
		plinkMarkerLoader.registerThread(thread);

		return plinkMarkerLoader;
	}



}
