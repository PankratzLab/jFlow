package one.JL;

import java.util.ArrayList;

import common.Array;
import common.Files;
import seq.manage.ReferenceGenome;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import filesys.Segment;

public class ReferenceGenomeTest {

	public static void testRef(Project proj) {

		// basic
		ReferenceGenome referenceGenome = new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), proj.getLog());
		String[] test1 = referenceGenome.getSequenceFor(new Segment((byte) 26, 1, 50));

		// other testing
		String[] test = referenceGenome.getSequenceFor(new Segment((byte) 26, 1, 50));
		String[] bah = Files.getFirstNLinesOfFile(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), 2, proj.getLog());
		System.out.println(Array.toStr(test));
		System.out.println(Array.toStr(bah));

		for (int i = 0; i < test.length; i++) {
			System.out.println(test[i] + "\t" + bah[1].charAt(i));// these should be equal
		}
		System.out.println(Array.toStr(test));
		System.out.println(Array.toStr(bah));
		ArrayList<Segment> tseg = new ArrayList<Segment>();

		for (int i = 1; i < 23; i++) {
			tseg.add(new Segment((byte) i, referenceGenome.getIndexedFastaSequenceFile().getSequenceDictionary().getSequence(i).getSequenceLength() - 10000, referenceGenome.getIndexedFastaSequenceFile().getSequenceDictionary().getSequence(i).getSequenceLength() - 9800));
		}

		int off = tseg.size();
		for (int i = 1; i < 23; i++) {
			Segment segment1 = tseg.get(i);
			Segment segment2 = tseg.get(off - i);

			tseg.add(segment1);
			tseg.add(segment2);
		}
		for (int i = 1; i < 23; i++) {
			tseg.add(new Segment((byte) i, 1, 100));
		}

		long time = System.currentTimeMillis();
		proj.getLog().reportTimeInfo("Starting da query");
		for (int i = 0; i < tseg.size(); i++) {
			long disTime = System.currentTimeMillis();
			System.out.println("On: " + tseg.get(i).getUCSClocation());
			String[] bases = referenceGenome.getSequenceFor(tseg.get(i));
			if (Array.unique(bases).length > 2) {
				System.out.println(tseg.get(i).getUCSClocation() + "\t" + Array.toStr(bases, ""));
			}
			proj.getLog().reportTimeElapsed(disTime);
		}
		proj.getLog().reportTimeElapsed(time);

		Segment[] markerSegs = new Segment[proj.getMarkerNames().length];
		MarkerSet markerSet = proj.getMarkerSet();
		for (int i = 0; i < markerSegs.length; i++) {
			markerSegs[i] = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i] - 50, markerSet.getPositions()[i] + 50);
		}
		proj.getLog().reportTimeInfo("Going for every marker in the project n=" + markerSegs.length);
		time = System.currentTimeMillis();

		for (int i = 0; i < markerSegs.length; i++) {
			if (i % 1000 == 0) {
				proj.getLog().reportTimeInfo(i + "");
			}
			referenceGenome.getSequenceFor(markerSegs[i]);
		}
		proj.getLog().reportTimeElapsed(time);

	}

	public static void testRef2(Project proj) {

		// basic
		ReferenceGenome referenceGenome = new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), proj.getLog());
		Segment[] markerSegs = new Segment[proj.getMarkerNames().length];
		MarkerSet markerSet = proj.getMarkerSet();
		for (int i = 0; i < markerSegs.length; i++) {
			markerSegs[i] = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i] - 50, markerSet.getPositions()[i] + 50);
		}
		proj.getLog().reportTimeInfo("Going for every marker in the project n=" + markerSegs.length);
		long time = System.currentTimeMillis();
		for (int i = 0; i < markerSegs.length; i++) {
			if (i % 1000 == 0) {
				proj.getLog().reportTimeInfo(i + "");
			}
			referenceGenome.getSequenceFor(markerSegs[i]);
		}
		proj.getLog().reportTimeElapsed(time);
		proj.getLog().reportTimeInfo(" n=" + markerSegs.length);
		
	}

	public static void main(String[] args) {
		Project proj = new Project(args[0], false);
		testRef2(proj);
	}
}
