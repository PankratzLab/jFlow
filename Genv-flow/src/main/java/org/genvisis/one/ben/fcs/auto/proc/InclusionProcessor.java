package org.genvisis.one.ben.fcs.auto.proc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class InclusionProcessor extends AbstractSampleProcessor {

	String outDir;

	public InclusionProcessor(String o) {
		outDir = o;
	}

	@Override
	public void processSample(SampleNode sn, Logger log) throws IOException {
		if (!Files.exists(sn.fcsFile)) {
			return;
		}
		loadPopsAndGates(sn);
		loadData(sn);

		long t1 = System.nanoTime();
		HashSet<Gate> allGates = new HashSet<>(sn.gating.gateMap.values());
		HashMap<String, boolean[]> incl = new HashMap<>();

		for (Gate g : allGates) {
			boolean[] gating = g.gate(d);
			incl.put(g.getName(), gating);
		}

		String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
		String outFile = outDir + cleanedName + "/" + cleanedName + ".incl.xln";
		new File(outDir + cleanedName + "/").mkdirs();

		PrintWriter writer = Files.getAppropriateWriter(outFile);
		writer.print("Index");
		for (Gate g : allGates) {
			writer.print("\t" + g.getName());
		}
		writer.println();

		for (int i = 0; i < d.getCount(); i++) {
			writer.print(Integer.toString(i));
			for (Gate g : allGates) {
				writer.print("\t" + Boolean.toString(incl.get(g.getName())[i]));
			}
			writer.println();
		}
		writer.close();
		log.reportTime("Processed " + sn.fcsFile + " in " + ext.getTimeElapsedNanos(t1));
	}
}
