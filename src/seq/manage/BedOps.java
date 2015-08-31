package seq.manage;

import filesys.Segment;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.IndexFactory.IndexType;

import java.io.File;
import java.io.IOException;

import common.Files;
import common.Logger;

/**
 * @author lane0212 <br>
 *         Class to handle some common .bed file ops
 */
public class BedOps {

	/**
	 * @param bedFile
	 * @param log
	 * @return whether the index file exists, or was able to be created
	 */
	public static boolean verifyBedIndex(String bedFile, Logger log) {
		boolean created = false;
		if (Files.exists(bedFile)) {
			String indexFile = bedFile + Tribble.STANDARD_INDEX_EXTENSION;
			if (Files.exists(indexFile)) {
				log.reportTimeInfo("Detected index file " + indexFile);
				created = true;
			} else {
				Index index = IndexFactory.createIndex(new File(bedFile), new BEDCodec(), IndexType.LINEAR);
				try {
					index.writeBasedOnFeatureFile(new File(bedFile));
				} catch (IOException e) {
					e.printStackTrace();
				}
				if (Files.exists(indexFile)) {
					created = true;
				}
			}
		} else {
			log.reportFileNotFound(bedFile);
		}
		return created;

	}

	public static Segment getSegment(BEDFeature bedFeature, Logger log) {
		return new Segment(bedFeature.getChr(), bedFeature.getStart(), bedFeature.getEnd());
	}

}
