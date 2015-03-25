package seq.manage;

import java.util.ArrayList;
import java.util.List;

import common.Logger;
import common.ext;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212
 *
 *         Flip a vcf's calls from alt to ref
 */
public class VCFFlip {
	// VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);

	private static VCFFileReader reader;

	public static void flip(String inputVCF, Logger log) {
		reader = new VCFFileReader(inputVCF, true);
		String output = inputVCF;
		if (inputVCF.endsWith(".vcf.gz")) {
			output = output.replaceAll(".vcf.gz", ".flipped.vcf.gz");
		} else {
			output = ext.addToRoot(inputVCF, ".flipped");
		}
		VariantContextWriter writer = VCFOps.initWriter(output, null, reader.getFileHeader().getSequenceDictionary());

		VariantContextWriterBuilder builderError = new VariantContextWriterBuilder().setOutputFile(inputVCF + ".errors.vcf.gz");
		builderError.setBuffer(2);
		builderError.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());
		VariantContextWriter writerError = builderError.build();

		VCFHeader vcfHeader = new VCFHeader(reader.getFileHeader());
		VCFFormatHeaderLine vcfFormatHeaderLine = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Ref to alt flipped genotypes from " + inputVCF);
		vcfHeader.addMetaDataLine(vcfFormatHeaderLine);
		writer.writeHeader(vcfHeader);
		writerError.writeHeader(vcfHeader);

		int skipped = 0;
		int count = 0;
		for (VariantContext vc : reader) {
			count++;
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			List<Allele> alleles = vc.getAlleles();
			List<Allele> newAlleles = new ArrayList<Allele>();
			Allele newRef = null;
			Allele newAlt = null;

			for (Allele a : alleles) {
				if (a.isReference()) {
					if (newAlt != null) {
						log.reportTimeError("Multiple new alts....");
						return;
					}
					newAlt = Allele.create(a.getBases(), false);
					newAlleles.add(Allele.create(a.getBases(), false));
				} else {
					if (newRef != null) {
						log.reportTimeError("Multiple new refs...." + newRef.getBaseString());
					}
					newRef = Allele.create(a.getBases(), true);

					newAlleles.add(Allele.create(a.getBases(), true));
				}
			}
			if (newRef.getBases().length != newAlt.getBases().length) {
				builder.start(vc.getStart());
				builder.stop(vc.getStart() + newRef.getBases().length - 1);

			}
			builder.alleles(newAlleles);
			ArrayList<Allele> newHomRef = new ArrayList<Allele>();
			newHomRef.add(newRef);
			newHomRef.add(newRef);
			ArrayList<Allele> newHomAlt = new ArrayList<Allele>();
			newHomAlt.add(newAlt);
			newHomAlt.add(newAlt);

			ArrayList<Allele> newHet = new ArrayList<Allele>();
			newHet.add(newRef);
			newHet.add(newAlt);

			GenotypesContext genotypesContext = vc.getGenotypes();
			GenotypesContext newGenotypes = GenotypesContext.create();
			for (Genotype g : genotypesContext) {
				GenotypeBuilder gbuilder = new GenotypeBuilder(g);
				if (g.isHomRef()) {
					gbuilder.alleles(newHomAlt);
				} else if (g.isHomVar()) {
					gbuilder.alleles(newHomRef);
				} else if (g.isHet()) {
					gbuilder.alleles(newHet);
				}
				newGenotypes.add(gbuilder.make());
			}
			builder.genotypes(newGenotypes);
			try {
				VariantContext flippedVC = builder.make();
				writer.add(flippedVC);

			} catch (IllegalArgumentException ile) {
				log.reportException(ile);
				skipped++;
				log.reportTimeError("Could not flip variant context " + vc.toStringWithoutGenotypes() + ", reverting to original context...");
				writerError.add(vc);
				writer.add(vc);
			}
			if (count % 100000 == 0) {
				log.reportTimeInfo("Flipped " + count + " variants");
			}

		}
		writer.close();
		writerError.close();
		if (skipped > 0) {
			log.reportTimeWarning("Some variants were invalid, we skipped " + skipped);
		}
		reader.close();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "VCFFlip.dat";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.manage.VCFFlip requires 0-1 arguments\n" + "   (1) vcf file (i.e. vcf=" + vcf + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			flip(vcf, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
