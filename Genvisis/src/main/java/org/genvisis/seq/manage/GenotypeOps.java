package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.Logger;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;

public class GenotypeOps {

	public static List<Allele> getNoCall() {
		ArrayList<Allele> noCall = new ArrayList<Allele>();
		noCall.add(Allele.NO_CALL);
		noCall.add(Allele.NO_CALL);
		return noCall;
	}

	public static String[] getGenoAnnotationsFor(String[] annosToGet, Genotype g, String defaultValue) {
		String[] annos = new String[annosToGet.length];
		for (int i = 0; i < annos.length; i++) {
			String tmp = defaultValue;
			if (g.hasAnyAttribute(annosToGet[i])) {
				tmp = g.getAnyAttribute(annosToGet[i]).toString();
			}
			annos[i] = tmp;
		}
		return annos;
	}

	public static String[][] getGenoFormatKeys(String vcf, Logger log) {
		VCFFileReader reader = new VCFFileReader(new File(vcf), false);
		String[] annotationKeys = new String[reader.getFileHeader().getFormatHeaderLines().size()];
		String[] descriptions = new String[reader.getFileHeader().getFormatHeaderLines().size()];
		int index = 0;
		for (VCFFormatHeaderLine vcfHeaderLine : reader.getFileHeader().getFormatHeaderLines()) {
			annotationKeys[index] = vcfHeaderLine.getID();
			descriptions[index] = vcfHeaderLine.getDescription();
			index++;
		}
		reader.close();
		return new String[][] { annotationKeys, descriptions };

	}

}
