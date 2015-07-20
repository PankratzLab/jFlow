package cnv.annotation;

import java.util.ArrayList;

import common.Positions;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import filesys.Segment;

public class LocusAnnotation {
	private String locusName;
	private Segment seg;
	private String ref;
	private String[] alts;

	private AnnotationData[] annotations;

	public VariantContext annotateLocus() {
		ArrayList<Allele> alleles = new ArrayList<Allele>();
		Allele refA = Allele.create(ref, true);
		alleles.add(refA);
		for (int i = 0; i < alts.length; i++) {
			Allele altA = Allele.create(alts[i], false);
			alleles.add(altA);
		}

		VariantContextBuilder vBuilder = new VariantContextBuilder();
		vBuilder.chr("chr"+seg.getChr());
		vBuilder.start(seg.getStart());
		vBuilder.stop(seg.getStop());
		vBuilder.id(locusName);
		vBuilder.alleles(alleles);
		if (annotations != null) {
			for (int i = 0; i < annotations.length; i++) {
				annotations[i].addAnnotation(vBuilder);
			}
		}
		return vBuilder.make();
	}

	public Segment getSeg() {
		return seg;
	}

	public AnnotationData[] getAnnotations() {
		return annotations;
	}

	public static class Builder {
		private String ref = "A";
		private String[] alts = new String[] { "N" };
		private AnnotationData[] annotations = null;

		public Builder ref(String ref) {
			this.ref = ref;
			return this;
		}

		public Builder alts(String[] alts) {
			this.alts = alts;
			return this;
		}

		public Builder annotations(AnnotationData[] annotations) {
			this.annotations = annotations;
			return this;
		}

		public LocusAnnotation build(String locusName, Segment seg) {
			return new LocusAnnotation(this, locusName, seg);
		}
	}

	private LocusAnnotation(Builder builder, String locusName, Segment seg) {
		this.locusName = locusName;
		this.seg = seg;
		this.annotations = builder.annotations;
		this.ref = builder.ref;
		this.alts = builder.alts;
	}
}
