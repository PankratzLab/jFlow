package org.genvisis.one.ben.flowannot;

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Annotator implements IAnnotator {

	Logger log = new Logger();

	ArrayList<String> fcsKeys = new ArrayList<>();
	HashMap<String, HashMap<String, AnnotatedImage>> imageMap = new HashMap<>();
	ArrayList<Annotation> annotations = new ArrayList<>();

	public HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap() {
		return imageMap;
	}

	@Override
	public void saveAnnotations(String annotFile) {
		// TODO Auto-generated method stub

	}

	@Override
	public void replaceAnnotation(Annotation prevAnnot, Annotation newAnnot) {
		// TODO Auto-generated method stub

	}

	@Override
	public void removeAnnotation(Annotation annotation) {
		// TODO Auto-generated method stub

	}

	@Override
	public void loadImgDir(String dir) {
		File dFil = new File(dir);
		File[] subDirs = dFil.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});
		for (File d : subDirs) {
			String fcsFilename = d.getName();
			fcsKeys.add(fcsFilename);
			HashMap<String, AnnotatedImage> fcsImgs = new HashMap<>();
			imageMap.put(fcsFilename, fcsImgs);
			String[] imgFiles = d.list(new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
					return name.startsWith(fcsFilename) && (name.endsWith(".png") || name.endsWith(".jpg"));
				}
			});
			for (String img : imgFiles) {
				int gateInd = -1;
				for (int i = 0; i < GateTree.GATE_DIMS.length; i++) {
					boolean containsAll = true;
					for (String s : GateTree.GATE_DIMS[i]) {
						if (!img.contains(s)) {
							containsAll = false;
							break;
						}
					}
					if (containsAll) {
						gateInd = i;
						break;
					}
				}
				if (gateInd >= 0) {
					AnnotatedImage ai = new AnnotatedImage(gateInd + "", gateInd == 0);
					ai.imageFile = ext.verifyDirFormat(d.getAbsolutePath()) + img;
					String name = ArrayUtils.toStr(GateTree.GATE_DIMS[gateInd], " v ");
					ai.gateName = name;
					fcsImgs.put(name, ai);
				}
			}
		}
	}

	@Override
	public void loadAnnotations(String annotFile) {
		// TODO Auto-generated method stub

	}

	@Override
	public String getLoadedAnnotationFile() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void deleteAnnotation(Annotation annotation) {
		// TODO Auto-generated method stub

	}

	@Override
	public void applyAnnotation(Annotation annotation) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addNewAnnotation(Annotation newAnnotation) {
		// TODO Auto-generated method stub

	}


}
