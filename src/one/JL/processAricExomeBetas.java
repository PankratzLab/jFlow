package one.JL;

import java.io.File;
import java.util.ArrayList;

import one.JL.BetaOptimizer.MarkerRsFormat;
import common.Array;
import common.Files;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.Resources.GENOME_BUILD;
import cnv.manage.Resources.GENOME_RESOURCE_TYPE;

public class processAricExomeBetas {

	public static void main(String[] args) {

		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		String out = proj.PROJECT_DIRECTORY.getValue() + "betaOpti/";
		new File(out).mkdirs();
		proj.AB_LOOKUP_FILENAME.setValue(out + "AB_LookupBeta.dat");
		ABLookup abLookup = new ABLookup();
		if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
			abLookup.parseFromAnnotationVCF(proj);
			abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(), proj.getLog());
		}
		MarkerSet markerSet = proj.getMarkerSet();
		abLookup = new ABLookup(markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog());
		String mapSer = out + "rsIDMap.ser";
		if (!Files.exists(mapSer)) {
			BetaOptimizer.mapToRsIds(markerSet, abLookup, GENOME_RESOURCE_TYPE.DB_SNP147.getResource(GENOME_BUILD.HG19).getResource(proj.getLog()), markerSet.getMarkerNames(), mapSer, proj.getLog());
		}
		ArrayList<MarkerRsFormat> markerRsFormats = MarkerRsFormat.readSerial(mapSer, proj.getLog());
		ArrayList<String> rsOut = new ArrayList<String>();
		rsOut.add("MarkerName\tPos\trsID\trsRef\trsAlt\tmarkerA\tmarkerB\tconfig\tsiteType");
		for (int i = 0; i < markerRsFormats.size(); i++) {
			MarkerRsFormat m = markerRsFormats.get(i);
			rsOut.add(m.getMarkerName() + "\t" + m.getPosMarker() + "\t" + m.getRs() + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getType());
		}
		Files.writeArrayList(rsOut, out + "rsmatch.txt");

	}

}
