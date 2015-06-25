package seq.manage;

import java.util.ArrayList;
import java.util.List;

import htsjdk.variant.variantcontext.Allele;

public class GenotypeOps {

	
	
	public static List<Allele> getNoCall() {
		ArrayList<Allele> noCall = new ArrayList<Allele>();
		noCall.add(Allele.NO_CALL);
		noCall.add(Allele.NO_CALL);
		return noCall;
	}

}
