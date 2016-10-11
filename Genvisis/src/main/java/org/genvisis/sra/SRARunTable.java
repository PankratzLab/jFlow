package org.genvisis.sra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;

/**
 * Class for managing SRA runtables - which can be used to obtain extra information from samples
 * downloaded from the short read archive
 *
 */
public class SRARunTable extends HashMap<String, SRASample> {

	/**
	 *
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Corresponds to the SRA sample ID
	 *
	 */
	private static final String RUN_S = "Run_s";

	/**
	 * The study ID upon submission to the SRA, useful to matching up to other data
	 */
	private static final String SUBMITTED_SAMPLE_ID_S = "submitted_sample_id_s";

	private static final String ASSAY_TYPE_S = "Assay_Type_s";

	private static final String ASSEMBLY_TYPE = "AssemblyName_s";
	private static final String PLATFORM_S = "Platform_s";

	/**
	 * @param log
	 */
	public SRARunTable() {
		super();
	}

	public String[] getAllRunSFiles() {
		ArrayList<String> sraFiles = new ArrayList<String>();
		for (String sample : keySet()) {
			sraFiles.add(get(sample).getRunS());
		}
		return Array.toStringArray(sraFiles);
	}

	/**
	 * @param sraTable load this run table
	 * @param log
	 * @return {@link SRARunTable} filled with {@link SRASample}
	 */
	public static SRARunTable load(String sraTable, Logger log) {
		SRARunTable sraRunTable = new SRARunTable();

		try {
			String[] requiredHeader = new String[] {RUN_S, SUBMITTED_SAMPLE_ID_S, ASSAY_TYPE_S,
																							ASSEMBLY_TYPE, PLATFORM_S};

			BufferedReader reader = Files.getAppropriateReader(sraTable);

			int[] indices = ext.indexFactors(	requiredHeader, reader.readLine().trim().split("\t"), true,
																				false);
			int numLoaded = 0;
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String[] parsed = Array.subArray(line, indices);

				String runS = parsed[0];
				String submittedSampleID = parsed[1];
				ASSAY_TYPE aType = ASSAY_TYPE.valueOf(parsed[2]);
				ASSEMBLY_NAME name = ASSEMBLY_NAME.valueOf(parsed[3].toUpperCase());
				PLATFORM platform = PLATFORM.valueOf(parsed[4]);

				if (platform == PLATFORM.ILLUMINA) {
					sraRunTable.put(runS, new SRASample(runS, submittedSampleID, name, aType, platform));
					numLoaded++;
				}
			}
			log.reportTimeInfo("Loaded " + numLoaded + " ILLUMINA samples from " + sraTable);
			reader.close();
		} catch (FileNotFoundException e) {
			log.reportException(e);
		} catch (IOException e) {
			log.reportException(e);

		}

		return sraRunTable;
	}

	public static void main(String[] args) {}

}
