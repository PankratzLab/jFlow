package org.genvisis.one.ben.flowannot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class Annotator implements IAnnotator {

	private ArrayList<String> fcsKeys = new ArrayList<>();
	private HashMap<String, HashMap<String, AnnotatedImage>> imageMap = new HashMap<>();
	private ArrayList<AnnotatedImage.Annotation> annotations = new ArrayList<>();

	// TODO add to this:
	private HashMap<AnnotatedImage.Annotation, ArrayList<AnnotatedImage>> annotMap = new HashMap<>();

	public ArrayList<String> getFCSKeys() {
		return fcsKeys;
	}

	public HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap() {
		return imageMap;
	}

	@Override
	public ArrayList<AnnotatedImage.Annotation> getAnnotations() {
		return this.annotations;
	}

	@Override
	public void addNewAnnotation(AnnotatedImage.Annotation newAnnotation) {
		this.annotations.add(newAnnotation);
	}

	@Override
	public void replaceAnnotation(AnnotatedImage.Annotation prevAnnot,
																AnnotatedImage.Annotation newAnnot) {
		this.annotations.set(this.annotations.indexOf(prevAnnot), newAnnot);
		for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
			for (AnnotatedImage ai : annMap.values()) {
				if (ai.getAnnotations().contains(prevAnnot)) {
					ai.getAnnotations().set(ai.getAnnotations().indexOf(prevAnnot), newAnnot);
				}
			}
		}
	}

	@Override
	public void deleteAnnotation(AnnotatedImage.Annotation annotation) {
		this.annotations.remove(annotation);
		for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
			for (AnnotatedImage ai : annMap.values()) {
				ai.getAnnotations().remove(annotation);
			}
		}
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
				String name = img.substring(fcsFilename.length() + 1, img.length() - 4);
				for (int i = 0; i < GateTree.GATE_TREE.length; i++) {
					if (GateTree.GATE_TREE[i][0].equals(name)
							|| ext.replaceWithLinuxSafeCharacters(GateTree.GATE_TREE[i][0]).equals(name)) {
						gateInd = i;
						break;
					}
				}
				if (gateInd >= 0) {
					AnnotatedImage ai = new AnnotatedImage(gateInd + "", gateInd == 0);
					ai.setImageFile(ext.verifyDirFormat(d.getAbsolutePath()) + img);
					// name = ArrayUtils.toStr(GateTree.GATE_DIMS[gateInd], " v ");
					ai.setGateName(GateTree.GATE_TREE[gateInd][0]);
					fcsImgs.put(GateTree.GATE_TREE[gateInd][0], ai);
				}
			}
		}
	}

	private static final String ANNOT_TOKEN = "@ANNOT";
	private static final String IMAGE_TOKEN = "@IMAGE";

	@Override
	public void loadAnnotations(String annotFile) throws IOException {
		BufferedReader reader = Files.getAppropriateReader(annotFile);
		String line = null;

		while ((line = reader.readLine()) != null) {
			if ("".equals(line))
				continue;
			if (line.startsWith(ANNOT_TOKEN)) {
				String[] pts = line.split("\t");
				AnnotatedImage.Annotation a = new AnnotatedImage.Annotation(pts[1], pts[2]);
				this.annotations.add(a);
			} else if (line.startsWith(IMAGE_TOKEN)) {
				String[] pts = line.split("\t")[1].split("\\|", -1);
				AnnotatedImage ai = new AnnotatedImage(pts[0], GateTree.GATE_TREE[0][0].equals(pts[0]));
				String imgFile = pts[1].equals("") ? null : pts[1];
				ai.setImageFile(imgFile);
				ai.setMissing(imgFile == null || !Files.exists(imgFile));
				if (pts.length > 2) {
					for (int i = 2; i < pts.length; i++) {
						for (AnnotatedImage.Annotation a : this.annotations) {
							if (a.annotation.equals(pts[i])) {
								ai.getAnnotations().add(a);
							}
						}
					}
				}
				String file = ext.removeDirectoryInfo(imgFile);
				String fcsFile = file.substring(0, file.indexOf(".fcs.") + 4);
				if (!fcsKeys.contains(fcsFile)) {
					fcsKeys.add(fcsFile);
				}
				HashMap<String, AnnotatedImage> map = imageMap.get(fcsFile);
				if (map == null) {
					map = new HashMap<>();
					imageMap.put(fcsFile, map);
				}
				map.put(ai.getGateName(), ai);
			}
		}
		reader.close();
	}

	public void saveAnnotation(AnnotatedImage.Annotation annot, String file) {
		PrintWriter writer = Files.getAppropriateWriter(file);
		writer.println(ANNOT_TOKEN + "=" + annot.annotation);
		writer.println();
		for (String fcs : fcsKeys) {
			for (Entry<String, AnnotatedImage> ent : imageMap.get(fcs).entrySet()) {
				if (ent.getValue().getAnnotations().contains(annot)) {
					writer.println(fcs + "\t" + ent.getKey());
				}
			}
		}
		writer.flush();
		writer.close();
	}

	@Override
	public void saveAnnotations(String annotFile) {
		HashMap<AnnotatedImage.Annotation, ArrayList<String>> map = new HashMap<>();
		for (AnnotatedImage.Annotation a : annotations) {
			map.put(a, new ArrayList<>());
		}
		for (String f : fcsKeys) {
			for (Entry<String, AnnotatedImage> ent : imageMap.get(f).entrySet()) {
				for (AnnotatedImage.Annotation a : ent.getValue().getAnnotations()) {
					map.get(a).add(ent.getKey());
				}
			}
		}
		PrintWriter writer = Files.getAppropriateWriter(annotFile);
		for (AnnotatedImage.Annotation a : map.keySet()) {
			StringBuilder sb = new StringBuilder(ANNOT_TOKEN).append("\t").append(a.annotation)
																											 .append("\t").append(a.mnemonic);
			writer.println(sb.toString());
		}
		writer.println();
		for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
			for (Entry<String, AnnotatedImage> ent : annMap.entrySet()) {
				writer.println(IMAGE_TOKEN + "\t" + ent.getValue().exportToString());
			}
		}
		writer.flush();
		writer.close();
	}

}
