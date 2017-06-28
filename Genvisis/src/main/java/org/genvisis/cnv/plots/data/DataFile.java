package org.genvisis.cnv.plots.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

import org.genvisis.CLI.Arg;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class DataFile {
	String filename;
	String sysFilename;
	Logger log;

	String[] fullHeader;
	String[] loadedHdr;
	String[][] linkers;
	int[] linkIndices;
	String[] sortColumns;
	Arg[] sortTypes;
	boolean[] reqHdr;
	ArrayList<String[]> dataRows;
	ArrayList<DataListener> listeners;

	Thread loadThread;

	volatile boolean success = false;
	String failMsg = "";

	public DataFile(String file, Logger log, String[][] linkersWithAlts, boolean[] requiredHeaders) {
		if (requiredHeaders.length != linkersWithAlts.length) {
			throw new IllegalArgumentException(
																				 "PROGRAM ERROR: Linkers and required headers objects must be the same length!");
		}
		this.filename = file;
		this.sysFilename = (new File(file)).getAbsolutePath();
		this.linkers = linkersWithAlts;
		this.reqHdr = requiredHeaders;
		this.log = log;
		this.listeners = new ArrayList<>();
		if (!Files.exists(file)) {
			return;
		}
		loadAndIndexHeader();
	}

	private void loadAndIndexHeader() {
		this.fullHeader = Files.getHeaderOfFile(filename, log);
		linkIndices = ext.indexFactors(linkers, fullHeader, false, true, false, log,
																	 false);
	}

	public boolean hasRequiredData() {
		for (int i = 0; i < reqHdr.length; i++) {
			if (reqHdr[i] && linkIndices[i] == -1) {
				return false;
			}
		}
		return true;
	}

	public String[] getColumns() {
		return fullHeader;
	}

	public String[] getUnlinkedColumns() {
		String[] columns = new String[fullHeader.length
																	- (linkIndices.length - ArrayUtils.countIf(linkIndices, -1))];
		int ind = 0;
		for (int i = 0; i < fullHeader.length; i++) {
			if (ext.indexOfInt(i, linkIndices) == -1) {
				columns[ind++] = fullHeader[i];
			}
		}

		return columns;
	}

	public void pingWhenReady(DataListener rddy) {
		this.listeners.add(rddy);
	}

	private void pingReady() {
		for (DataListener dl : this.listeners) {
			dl.ping(this);
		}
	}

	public void loadData(HashMap<String, DataPipe> selData, boolean waitForData) {
		if (waitForData) {
			actuallyLoadData(selData, false);
		} else {
			loadThread = new Thread(() -> {
				actuallyLoadData(selData, true);
			}, "DataFile_" + filename);
			loadThread.setDaemon(true);
			loadThread.start();
		}
	}

	private HashMap<String, Integer> colIndices = new HashMap<>();

	public int getIndexOfData(String col) {
		if (col == null || "".equals(col) || !colIndices.containsKey(col)) {
			return -1;
		}
		return colIndices.get(col);
	}

	public void setSortColumns(String[] columnNames, Arg[] args) {
		this.sortColumns = columnNames;
		this.sortTypes = args;
	}

	private void actuallyLoadData(HashMap<String, DataPipe> selData, boolean ping) {
		log.reportTime("Loading file " + filename);
		this.dataRows = new ArrayList<>();
		int[] columnsToLoad = new int[selData.size()];
		int[] sortIndices = new int[sortColumns == null ? 0 : sortColumns.length];
		loadedHdr = new String[selData.size()];
		DataPipe[] pipes = new DataPipe[selData.size()];
		int ind = 0;
		for (int i = 0; i < fullHeader.length; i++) {
			if (selData.containsKey(fullHeader[i])) {
				columnsToLoad[ind] = i;
				colIndices.put(fullHeader[i], ind);
				pipes[ind] = selData.get(fullHeader[i]);
				loadedHdr[ind] = fullHeader[i];
				ind++;
			}
		}
		for (int i = 0; i < sortIndices.length; i++) {
			sortIndices[i] = ext.indexOfStr(sortColumns[i], loadedHdr);
			if (sortIndices[i] == -1) {
				log.reportTime("ERROR - Sort column [" + sortColumns[i]
											 + "] was not included in loaded data or doesn't exist in the file.");
				log.report("Loaded Columns: " + ArrayUtils.toStr(loadedHdr, ", "));
				log.report("All Columns: " + ArrayUtils.toStr(fullHeader, ", "));
			}
		}

		int nLines = Files.countLines(filename, 1);
		log.reportTime("Reading " + nLines + " lines of data from file: " + filename);
		try {
			BufferedReader reader = Files.getAppropriateReader(filename);
			String line = reader.readLine(); // skip header
			String delim = ext.determineDelimiter(line.trim());
			String[] parts = null;
			while ((line = reader.readLine()) != null) {
				parts = line.trim().split(delim, -1);
				String[] data = new String[columnsToLoad.length];
				boolean validLine = true;
				for (int i = 0; i < data.length; i++) {
					String raw = parts[columnsToLoad[i]];
					data[i] = pipes[i] == null ? raw : pipes[i].pipe(raw);
					if (data[i] == null) {
						validLine = false;
						break;
					}
				}
				if (validLine) {
					this.dataRows.add(data);
				}
			}
			reader.close();
			log.reportTime("Loaded " + dataRows.size() + " lines of data from file: " + filename);
			if (sortIndices != null && sortIndices.length > 0) {
				this.dataRows.sort(new Comparator<String[]>() {
					public int compare(String[] o1, String[] o2) {
						for (int i = 0; i < sortIndices.length; i++) {
							int ind = sortIndices[i];
							int res = 0;
							if (sortTypes != null) {
								switch (sortTypes[i]) {
									case NUMBER:
										res = Double.compare(Double.parseDouble(o1[ind]), Double.parseDouble(o2[ind]));
										break;
									case STRING:
										res = o1[ind].compareTo(o2[ind]);
										break;
									case FILE:
									default:
										break;
								}
							}
							if (res != 0) {
								return res;
							}
						}
						return 0;
					};
				});
				log.reportTime("Sorted data from file: " + filename);
			}
			success = true;
		} catch (IOException e) {
			success = false;
			failMsg = e.getMessage();
		}
		if (ping) {
			pingReady();
		}
	}

	public int[] getLinkIndices() {
		return linkIndices;
	}

	public ArrayList<String[]> getAllData() {
		return dataRows;
	}

	public String getLinkedColumnName(int linkIndex) {
		int ind = getLinkIndices()[linkIndex];
		if (ind == -1) {
			return null;
		}
		return getColumns()[ind];
	}

	public boolean hasLinkedColumn(int linker) {
		return getLinkIndices()[linker] != -1;
	}

	public String getFilename() {
		return filename;
	}

	public boolean isLoaded() {
		return success;
	}

}
