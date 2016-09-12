package org.genvisis.one.JL.aging;

import java.io.IOException;
import java.util.ArrayList;

import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BamOps;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Trouble shooting some issues with build 37, etc.... chr XY and chr0 causing problems
 *
 */
public class TestSeqDictionary {

	public static void main(String[] args) {
		String bam = "/Volumes/Beta/data/aric_sra/test/bams/SRR1700989.bam";
		ArrayList<Segment> testSegs = new ArrayList<Segment>();
		for (int i = 0; i < 27; i++) {
			testSegs.add(new Segment((byte) i, 1, 1000));

		}
		SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);

		for (Segment cs : testSegs) {
			String chr = Positions.getChromosomeUCSC(cs.getChr(), false, true);
			if (cs.getChr() == 26) {
				chr = "MT";
			}
			new Logger().reportTime(chr);
			CloseableIterator<SAMRecord> iterator = reader.queryOverlapping(chr, cs.getStart(), cs.getStop());
			int start = 0;
			while (iterator.hasNext()) {
				SAMRecord record = iterator.next();
				if (start == 0) {
					new Logger().reportTimeInfo(record.toString());
				}
				start++;

			}
			iterator.close();
		}
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
