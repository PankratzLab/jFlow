package org.genvisis.seq.manage.vcfPile;

import java.io.File;
import java.io.PrintWriter;
import java.util.Iterator;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 Summarizes variants by region in a vcf
 */
public class VCFPile<T extends Segment> implements Iterator<PiledVcfRegion<T>> {
	private final String vcfFile;
	private final String[] samplesTopile;
	private final ReferenceGenome referenceGenome;
	private final RegionIteratorVCF<T> rIteratorVCF;
	private final Logger log;

	public VCFPile(String vcfFile, ReferenceGenome referenceGenome, String[] samplesTopile,
								 LocusSet<T> toPile, Logger log) {
		super();
		this.vcfFile = vcfFile;
		this.samplesTopile = samplesTopile;
		this.rIteratorVCF = new RegionIteratorVCF<T>(vcfFile, toPile, log);
		this.referenceGenome = referenceGenome;
		this.log = log;
	}

	public String getVcfFile() {
		return vcfFile;
	}

	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
		return rIteratorVCF.hasNext();
	}

	@Override
	public PiledVcfRegion<T> next() {
		PiledVcfRegion<T> pRegion = new PiledVcfRegion<T>(rIteratorVCF.getCurrentIndex(),
																											referenceGenome, samplesTopile);
		VariantContext[] vcs = rIteratorVCF.next();
		for (VariantContext vc : vcs) {
			pRegion.addVariantContext(vc, log);
		}
		return pRegion;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
	}

	public static void pileVCF(String vcfFile, String referenceGenomeFile, String regionsFile,
														 String outputDirectory) {
		if (outputDirectory == null) {
			outputDirectory = ext.parseDirectoryOfFile(vcfFile);
		}
		new File(outputDirectory).mkdirs();
		Logger log = new Logger(outputDirectory + VCFOps.getAppropriateRoot(vcfFile, true)
														+ ".vcfPile.log");

		log.reportTimeInfo("Loading regions from " + regionsFile);
		Segment[] segs = Segment.loadRegions(regionsFile, 0, 1, 2, false);
		LocusSet<Segment> toPile = new LocusSet<Segment>(segs, true, log) {

			/**
			 *
			 */
			private static final long serialVersionUID = 1L;
		};
		log.reportTimeInfo("Loading " + segs.length + " regions from ");

		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFile, log);
		String[] samples = VCFOps.getSamplesInFile(new VCFFileReader(new File(vcfFile), true));
		VCFPile<Segment> vcfPile = new VCFPile<Segment>(vcfFile, referenceGenome, samples, toPile, log);
		String output = outputDirectory + VCFOps.getAppropriateRoot(vcfFile, true)
										+ ".vcfPile.summary.txt";
		try {

			PrintWriter writer = Files.openAppropriateWriter(output);
			writer.println("REGION\tAVG_NUM_VAR\tAVG_GC\tAVG_DP\tAVG_GQ");
			int index = 0;
			while (vcfPile.hasNext()) {
				index++;
				PiledVcfRegion<Segment> pRegion = vcfPile.next();
				log.reportTimeInfo("On Region " + index + "\t" + pRegion.getRegion().getUCSClocation());
				String out = "";
				out += pRegion.getRegion().getUCSClocation();
				out += "\t" + ArrayUtils.mean(pRegion.getTotalCalledVar());
				out += "\t" + pRegion.getAvgGC();
				out += "\t" + ArrayUtils.mean(pRegion.getAvgDP());
				out += "\t" + ArrayUtils.mean(pRegion.getAvgGQ());
				writer.println(out);
				writer.flush();
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcfFile = "aVcf.vcf.gz";
		String referenceGenomeFile = "aref.fa";
		String regionsFile = "AgilentCaptureRegions.bed";

		String outputDirectory = null;

		// String logfile = null;
		// Logger log;

		String usage = "\n" + "seq.manage.vcfPile.VCFPile requires 0-1 arguments\n";
		usage += "   (1) vcf  (i.e. vcf=" + vcfFile + " (default))\n" + "";
		usage += "   (2) reference genome  (i.e. ref=" + referenceGenomeFile + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, "");
		usage += "   (4) a regions file (i.e. regions=" + regionsFile + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcfFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("ref=")) {
				referenceGenomeFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("regions=")) {
				regionsFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDirectory = arg.split("=")[1];
				numArgs--;
			}
			// else if (args[i].startsWith("log=")) {
			// logfile = args[i].split("=")[1];
			// numArgs--;
			// }
			else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// log = new Logger(logfile);
			pileVCF(vcfFile, referenceGenomeFile, regionsFile, outputDirectory);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
