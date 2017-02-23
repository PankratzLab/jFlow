package org.genvisis.one.ben;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeader.HEADER_FIELDS;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.seq.manage.ReferenceGenome;

import ca.mcgill.mcb.pcingola.genotypes.Genotypes;

public class VCFConstructor {
	
	Logger log;
	String inputFile;
	String dosageFile;
	String outputFile;
	
	String[] markers;
	int[][] locations;
	Map<String, List<Allele>> alleles;
	Map<String, Collection<Genotype>> genotypes;
	String[] ids;
	
	private void readInputFile() {
		// TODO check inputFile is set and exists and is readable
		markers = HashVec.loadFileToStringArray(inputFile, false, false, new int[] {0}, true, false, "\t");
		SnpMarkerSet markerSet = new SnpMarkerSet(markers);
		markerSet.parseSNPlocations(log);
		locations = markerSet.getChrAndPositionsAsInts();
		parseAlleles(markerSet.getMarkerNames());
	}
	
	private void readDosageFile() {
		// TODO check dosageFile is set and exists and is readable
		String[][] matr = HashVec.loadFileToStringMatrix(dosageFile, false, null, false);
		ids = ArrayUtils.subArray(Matrix.extractColumn(matr, 1), 1);
		String[] snps = ArrayUtils.subArray(matr[0], 2);
		
		genotypes = new HashMap<String, Collection<Genotype>>();
		for (int i = 0; i < snps.length; i++) {
			String mkr = snps[i];
			List<Allele> all = alleles.get(mkr);
			int m = i + 2;
			Collection<Genotype> genos = new ArrayList<Genotype>();
			for (int id = 0; id < ids.length; id++) {
				int ind = id + 1;
				String indDose = matr[ind][m];
				Genotype g;
				switch(indDose.charAt(0)) {
					case '0':
						g = GenotypeBuilder.create(ids[id], ArrayUtils.toList(new Allele[]{all.get(0), all.get(0)}));
						break;
					case '1':
						g = GenotypeBuilder.create(ids[id], ArrayUtils.toList(new Allele[]{all.get(0), all.get(1)}));
						break;
					case '2':
						g = GenotypeBuilder.create(ids[id], ArrayUtils.toList(new Allele[]{all.get(1), all.get(1)}));
						break;
					default:
						g = GenotypeBuilder.create(ids[id], new ArrayList<Allele>());
						break;
				}
				genos.add(g);
			}
			genotypes.put(mkr, genos);
		}
		
	}
	
	private void parseAlleles(String[] markerNames) {
		alleles = new HashMap<String, List<Allele>>();
		for (int i = 0; i < markerNames.length; i++) {
			if (markerNames[i].startsWith("rs")) {
				log.reportError("Cannot parse alleles from RS marker");
				alleles.put(markerNames[i], new ArrayList<Allele>());
			} else {
				String[] pts = markerNames[i].split(":");
				ArrayList<Allele> a = new ArrayList<Allele>();
				a.add(Allele.create(pts[2], true));
				a.add(Allele.create(pts[3], false));
				alleles.put(markerNames[i], a);
			}
		}
	}
	
	public void build() {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputFile);
		builder.clearOptions();
		builder.setOption(Options.INDEX_ON_THE_FLY);
		VCFHeader vcfHeader = new VCFHeader();

		
		SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(Resources.genome(GENOME_BUILD.HG19, log).getFASTA().getAbsolute(), log)
				.getIndexedFastaSequenceFile().getSequenceDictionary();

		builder.setReferenceDictionary(samSequenceDictionary);
		vcfHeader.setSequenceDictionary(samSequenceDictionary);
		VariantContextWriter writer = builder.build();
		writer.writeHeader(vcfHeader);

		for (int i = 0; i < markers.length; i++) {/// your data
			VariantContextBuilder builderVc = new VariantContextBuilder();
			builderVc.chr("chr" + locations[i][0]);
			int len = alleles.get(markers[i]).get(0).length() - 1;
			builderVc.alleles(alleles.get(markers[i]));
			builderVc.start(locations[i][1]);
			builderVc.stop(locations[i][1] + len);
			builderVc.id(markers[i]);
			builderVc.genotypes(genotypes.get(markers[i]));

			writer.add(builderVc.make());
		}
		writer.close();
	}
	
	public static void main(String[] args) {
		VCFConstructor c = new VCFConstructor();
		c.inputFile = "F:/temp/variantviewer/inputData.txt";
		c.outputFile = "F:/temp/variantviewer/data.vcf";
		c.dosageFile = "F:/temp/variantviewer/FRZ.db.xln";
		c.log = new Logger();
		c.readInputFile();
		c.readDosageFile();
		c.build();
	}
	
}
