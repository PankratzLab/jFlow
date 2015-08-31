package seq.qc;

import java.util.ArrayList;

import cnv.var.CNVariant;
import cnv.var.LocusSet;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.bed.BEDFeature;
import seq.manage.BEDFileReader;
import seq.manage.BedOps;
import common.Logger;
import common.Positions;
import common.ext;

public class Mappability<SEGMENT extends Segment> {

	private LocusSet<SEGMENT> set;
	private ArrayList<MappabilityResult<SEGMENT>> mappabilityResults;
	private String mappabilityFile;
	private Logger log;

	public Mappability(LocusSet<SEGMENT> set, String mappabilityFile, Logger log) {
		super();
		this.set = set;
		this.mappabilityFile = mappabilityFile;
		this.mappabilityResults = null;
		this.log = log;
	}

	public void computeMappability() {
		if (BedOps.verifyBedIndex(mappabilityFile, log)) {
			this.mappabilityResults = new ArrayList<Mappability.MappabilityResult<SEGMENT>>();
			BEDFileReader reader = new BEDFileReader(mappabilityFile, true);
			for (int i = 0; i < set.getLoci().length; i++) {
				if (i % 100 == 0) {
					log.reportTimeInfo(i + " of " + set.getLoci().length);
				}
				mappabilityResults.add(new MappabilityResult<SEGMENT>(reader, set.getLoci()[i], log));
			}
			reader.close();
		}
	}

	public static void computeCNVMappability(String mappabilityFile, String cnvFile, String geneTrackFile, Logger log) {
		BedOps.verifyBedIndex(mappabilityFile, log);
		LocusSet<GeneData> gLocusSet = GeneTrack.load(geneTrackFile, false).convertToLocusSet(log);
		CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile, false);
		LocusSet<CNVariant> cLocusSet = new LocusSet<CNVariant>(cnvs, true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		Mappability<CNVariant> cnMappability = new Mappability<CNVariant>(cLocusSet, mappabilityFile, log);
		cnMappability.computeMappability();

//		Mappability<GeneData> gMappability = new Mappability<GeneData>(gLocusSet, mappabilityFile, log);
//		gMappability.computeMappability();

	}

	private static class MappabilityResult<SEGMENT extends Segment> {
		private Logger log;
		private double cumulativeMapScore;
		private int numBases;
		private double averageMapScore;
		private SEGMENT t;

		public MappabilityResult(BEDFileReader reader, SEGMENT t, Logger log) {
			super();
			this.t = t;
			this.log = log;
			computeMappability(reader);
		}

		private void computeMappability(BEDFileReader reader) {
			CloseableIterator<BEDFeature> iterator = reader.query(Positions.getChromosomeUCSC(t.getChr(), true), t.getStart(), t.getStop());
			this.numBases = 0;
			double currentScore = -1;
			while (iterator.hasNext()) {
				BEDFeature bedFeature = iterator.next();
				Segment bedSeg = BedOps.getSegment(bedFeature, log);
				double mapScore = Double.NaN;
				try {
					mapScore = Double.parseDouble(bedFeature.getName());

				} catch (NumberFormatException nfe) {
					String error = "Could not convert " + bedFeature.getName() + " to mappability score";
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
				if (mapScore < 0) {
					String error = "Could not convert " + bedFeature.getName() + " to mappability score, negative mapScore";
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
				Segment union = bedSeg.getUnion(t, log);
				numBases += union.getSize();
				mapScore *= union.getSize();
				if (currentScore > 0) {
					currentScore += mapScore;
				} else {
					currentScore = mapScore;
				}
			}
			if (numBases != t.getSize()) {
				String warning = "Did not find scores for all bases in " + t.toAnalysisString() + " found " + numBases + " and should have found " + t.getSize();
				warning += "Setting remaining bases to score of 0";
				log.reportTimeError(warning);
				numBases += t.getSize() - numBases;
			}
			this.cumulativeMapScore = currentScore;
			this.averageMapScore = (double) cumulativeMapScore / numBases;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String mappabilityFile = "Mappability.bed";
		String cnvFile = "cnvs.cnv";
		String geneTrackFile = "RefSeq_hg19.gtrack";
		String usage = "\n" + "one.JL.Mappability requires 0-1 arguments\n";
		usage += "   (1) mappability file (i.e. mapFile=" + mappabilityFile + " (default))\n" + "";
		usage += "   (2) cnv file (i.e. cnvs=" + cnvFile + " (default))\n" + "";
		usage += "   (3) geneTrackFile  (i.e. genes=" + cnvFile + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("mapFile=")) {
				mappabilityFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvFile=")) {
				cnvFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("genes=")) {
				geneTrackFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Logger log = new Logger(ext.rootOf(cnvFile, false) + ".mappability.log");

			computeCNVMappability(mappabilityFile, cnvFile, geneTrackFile, log);
		} catch (Exception e) {

			e.printStackTrace();
		}
	}

}
