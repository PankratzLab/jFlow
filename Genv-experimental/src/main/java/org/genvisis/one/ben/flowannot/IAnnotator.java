package org.genvisis.one.ben.flowannot;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.HashMap;

import javax.imageio.ImageIO;
import javax.swing.UIManager;

import org.genvisis.common.Files;
import org.genvisis.common.ext;


public interface IAnnotator {

	static class Annotation {
		final String annotation;

		public Annotation(String ann) {
			this.annotation = ann;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((annotation == null) ? 0 : annotation.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Annotation other = (Annotation) obj;
			if (annotation == null) {
				if (other.annotation != null)
					return false;
			} else if (!annotation.equals(other.annotation))
				return false;
			return true;
		}

	}

	static class AnnotatedImage {
		private String gateName;
		private String imageFile;
		private SoftReference<BufferedImage> image;
		private ArrayList<Annotation> annots;
		private final boolean isRoot;
		private boolean missing = false;

		public AnnotatedImage(String ident, boolean isRoot) {
			this.gateName = ident;
			this.isRoot = isRoot;
			annots = new ArrayList<>();
		}

		public void setMissing(boolean miss) {
			missing = miss;
		}

		@Override
		public String toString() {
			return (isRoot && imageFile != null ? ext.rootOf(imageFile) + " | " : "")
						 + gateName + (missing ? " - MISSING" : "");
		}

		public BufferedImage getImage() {
			if (image == null) {
				createImage();
			}
			BufferedImage bi;
			while ((bi = image.get()) == null) {
				createImage();
			}
			return bi;
		}

		private void createImage() {
			if (image == null) {
				if (imageFile != null) {
					try {
						image = new SoftReference<BufferedImage>(ImageIO.read(new File(imageFile)));
					} catch (IOException e) {
						e.printStackTrace();
						image = new SoftReference<BufferedImage>(createIOExceptionImage(e));
					}
				} else {
					if (imageFile == null) {
						image = new SoftReference<BufferedImage>(createNoFileImage());
					} else if (!Files.exists(imageFile)) {
						image = new SoftReference<BufferedImage>(createMissingFileImage(imageFile));
					}
				}
			}
		}

		private static BufferedImage createImage(String msg, String msg2) {
			BufferedImage bi = new BufferedImage(800, 600, BufferedImage.TYPE_INT_RGB);
			Graphics g = bi.getGraphics();
			g.setColor(UIManager.getColor("Panel.background"));
			g.fillRect(0, 0, bi.getWidth(), bi.getHeight());
			g.setColor(Color.BLACK);
			g.setFont(g.getFont().deriveFont(fontSize));
			FontMetrics fm = g.getFontMetrics();
			g.drawString(msg, (bi.getWidth() / 2) - fm.stringWidth(msg) / 2,
									 bi.getHeight() / 2 - fm.getHeight());
			if (msg2 != null) {
				g.drawString(msg2, (bi.getWidth() / 2) - fm.stringWidth(msg2) / 2,
										 bi.getHeight() / 2 + ((int) (fm.getHeight() * 1.5)));
			}
			return bi;
		}

		static float fontSize = 16f;

		private BufferedImage createIOExceptionImage(IOException e) {
			return createImage("Exception when loading image:", e.getMessage());
		}

		public static BufferedImage createReadyImage() {
			return createImage("Select a gate to begin.", null);
		}

		private BufferedImage createNoFileImage() {
			return createImage("No image file found!", null);
		}

		private BufferedImage createMissingFileImage(String file) {
			return createImage("File missing:", file);
		}


		public boolean isRoot() {
			return isRoot;
		}

		public String getGateName() {
			return gateName;
		}

		public void setGateName(String name) {
			this.gateName = name;
		}

		public void setImageFile(String image) {
			this.imageFile = image;
		}

		public ArrayList<Annotation> getAnnotations() {
			return this.annots;
		}

	}



	void loadImgDir(String dir);

	void loadAnnotations(String annotFile);

	void saveAnnotations(String annotFile);

	void addNewAnnotation(Annotation newAnnotation);

	void deleteAnnotation(Annotation annotation);

	void replaceAnnotation(Annotation prevAnnot, Annotation newAnnot);

	String getLoadedAnnotationFile();

	void removeAnnotation(Annotation annotation);

	HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap();

	ArrayList<String> getFCSKeys();

	ArrayList<Annotation> getAnnotations();

}
