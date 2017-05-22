package org.genvisis.one.ben.flowannot;

import java.awt.image.BufferedImage;
import java.util.ArrayList;

import org.genvisis.common.ext;


public interface IAnnotator {

	static class Annotation {
		String annotation;

		// TODO hashcode/equals
	}

	static class AnnotatedImage {
		String gateIdent;
		String gateName;
		String imageFile;
		BufferedImage image;
		ArrayList<Annotation> annots;
		final boolean isRoot;
		boolean missing = false;

		public AnnotatedImage(String ident, boolean isRoot) {
			this.gateIdent = ident;
			this.gateName = ident;
			this.isRoot = isRoot;
			annots = new ArrayList<>();
		}

		public void setMissing(boolean miss) {
			missing = miss;
		}

		@Override
		public String toString() {
			return (isRoot && gateIdent != null && !gateIdent.equals(gateName) ? ext.rootOf(imageFile)
																																					 + " | "
																																				: "")
						 + gateName + (missing ? " - MISSING" : "");
		}

		public BufferedImage getImage() {
			// TODO SoftReference for image

			// TODO Auto-generated method stub
			return null;
		}
	}



	void loadImgDir(String dir);

	void loadAnnotations(String annotFile);

	void saveAnnotations(String annotFile);

	void addNewAnnotation(Annotation newAnnotation);

	void deleteAnnotation(Annotation annotation);

	void replaceAnnotation(Annotation prevAnnot, Annotation newAnnot);

	String getLoadedAnnotationFile();

	void applyAnnotation(Annotation annotation);

	void removeAnnotation(Annotation annotation);

}
