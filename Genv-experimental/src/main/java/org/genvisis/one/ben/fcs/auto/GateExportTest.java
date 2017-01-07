package org.genvisis.one.ben.fcs.auto;

import java.io.IOException;

import org.genvisis.common.Logger;
import org.genvisis.one.ben.fcs.auto.LeafDataSampler.GateAssignmentFactory;

public class GateExportTest {

//	private static final String WSP = "F:/Flow/gatingTest/";
//	private static final String FCS = "F:/Flow/gatingTest/2016-07-13_PANEL 1_ZF_Group one_F1631196_002.fcs";
	
	private static final String FCS = "/home/pankrat2/shared/flow/firstSetOf102/";
	private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	
	
	private void run() throws IOException {
		
		ProcessorFactory<? extends SampleProcessor> pf;
		
		pf = new GateAssignmentFactory();
//		pf = new PercentageWriterFactory();
//		pf = new LeafDataSamplerFactory("/scratch.global/cole0482/FCS/", new Logger());
		
		SamplingPipeline sp = new SamplingPipeline(null, WSP, FCS, null, "/scratch.global/cole0482/FCS/", pf);
		
		sp.run();
		
		while (!sp.checkFinished()) {
			Thread.yield();
		}
		
	}

	public static void main(String[] args) throws IOException {
		new GateExportTest().run();
	}
	
}
