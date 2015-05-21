package seq.manage.vcfPile;

import filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.util.ArrayList;
import java.util.Iterator;

import seq.manage.VCOps;
import cnv.var.LocusSet;
import common.Logger;
import common.Positions;

/**
 * @author lane0212 Iterate by region over a vcf;
 * 
 */
public class RegionIteratorVCF<T extends Segment> implements Iterator<VariantContext[]> {
	private String vcfFile;// force the file name to prevent reader being passed across threads
	private VCFFileReader reader;
	private LocusSet<T> set;
	private int setIndex;
	private Logger log;

	public RegionIteratorVCF(String vcfFile, LocusSet<T> segs, Logger log) {
		super();
		this.reader = new VCFFileReader(vcfFile, true);
		this.setIndex = 0;
		this.set = segs;
		this.log = log;
	}

	@Override
	public boolean hasNext() {
		boolean hasNext = setIndex < set.getLoci().length;
		if (!hasNext) {
			reader.close();
		}
		return hasNext;
	}

	public T getCurrentIndex() {
		return set.getLoci()[setIndex];
	}

	@Override
	public VariantContext[] next() {
		T region = set.getLoci()[setIndex];
		setIndex++;
		CloseableIterator<VariantContext> rIter = reader.query(Positions.getChromosomeUCSC(region.getChr(), true), region.getStart(), region.getStop());
		ArrayList<VariantContext> tmp = new ArrayList<VariantContext>();
		while (rIter.hasNext()) {
			VariantContext vc = rIter.next();// apparently an approximate query
			if (VCOps.getSegment(vc).overlaps(region)) {
				tmp.add(vc);
			}
		}
		return tmp.toArray(new VariantContext[tmp.size()]);
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
	}

}
