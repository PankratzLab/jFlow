package org.genvisis.one.george;

public class PipelineLauncher {
	public static void main(String[] args) {
		GMMPipeline pipeline = new GMMPipeline("./fcs");
		pipeline.run();
	}
}
