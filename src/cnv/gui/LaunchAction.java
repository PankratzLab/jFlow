package cnv.gui;

import cnv.plots.ScatterPlot;
import cnv.plots.Trailer;
import common.Array;
import common.ext;
import common.Files;
import java.awt.event.ActionEvent;
import javax.swing.AbstractAction;
import java.awt.Color;
import java.awt.Toolkit;

import javax.swing.*;

import cnv.filesys.Project;
import cnv.filesys.Sample;

public class LaunchAction extends AbstractAction {
	public static final long serialVersionUID = 1L;
	public static final int LAUNCH_TRAILER = 1;
	public static final int COPY_ID = 2;
	public static final int APPEND_ID = 3;
	public static final int APPEND_ID_TO_FILE = 4;
	public static final int LAUNCH_SCATTER = 5;
	public static final int LABEL_ONLY = 6;

	private Project proj;
	private String sample;
	private String filename;
	private String marker;
	private String[] loc;
	private int type;
	private boolean jar;
	private int plotStartX;//zx
	private int[] plotStartY;//zx
	private int plotWidth;//zx
	private int plotHeight;//zx

	//added by Zack.
	public LaunchAction(Project proj, String sample, String[] loc, Color color) {
		super(sample+" "+Array.toStr(loc, " / "));
		this.type = LAUNCH_TRAILER;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		this.loc = loc;
		this.plotStartX = Trailer.DEFAULT_STARTX; //zx
		this.plotStartY = new int[loc.length]; //zx
		this.plotWidth = Toolkit.getDefaultToolkit().getScreenSize().width - 30 - Trailer.DEFAULT_STARTX;
		this.plotHeight = (Toolkit.getDefaultToolkit().getScreenSize().height-50)/loc.length; //zx
		for (int i=0; i<loc.length; i++) {
			this.plotStartY[i] = 1 + i*(Toolkit.getDefaultToolkit().getScreenSize().height-50)/loc.length; //zx
		}
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}
	
	public LaunchAction(Project proj, String sample, String loc, Color color) {
		super(sample+" "+loc);
		this.type = LAUNCH_TRAILER;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		this.loc = new String[] {loc};
		this.plotStartX = Trailer.DEFAULT_STARTX;
		this.plotStartY = new int[] {Trailer.DEFAULT_STARTY};//zx
		this.plotWidth = Toolkit.getDefaultToolkit().getScreenSize().width - 30 - Trailer.DEFAULT_STARTX;
		this.plotHeight = Trailer.DEFAULT_HEIGHT;//zx
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}
	
	public LaunchAction(Project proj, String marker, Color color) {
		super(marker);
		this.type = LAUNCH_SCATTER;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.marker = marker;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}
	
	public LaunchAction(int type, Project proj, String sample, Color color) {
		super(sample);
		this.type = type;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}

	public LaunchAction(String filename, Project proj, String sample, Color color) {
		super(sample);
		this.type = APPEND_ID_TO_FILE;
		this.filename = filename;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}

	public LaunchAction(String text) {
		super(text);
		this.type = LABEL_ONLY;
	}
	
	public LaunchAction(String text, boolean append) {
		super((append?"Append":"Copy")+" to clipboard: "+ext.replaceAllWith(text, "\t", "  "));
		if (append) {
			this.type = APPEND_ID;
		} else {
			this.type = COPY_ID;
		}
		sample = text;
	}

	public void actionPerformed(ActionEvent e) {
		switch (type) {
        case LAUNCH_TRAILER:
    		ext.setClipboard(sample+"\t"+loc);
    		for (int i = 0; i < loc.length; i++) {
        		new Trailer(proj,
        					sample,
        					proj.getFilenames(Project.CNV_FILENAMES),
        					loc[i].endsWith("p")||loc[i].endsWith("q")?loc[i].substring(0, loc[i].length()-1):loc[i],
        					plotStartX,
        					plotStartY[i],
        					plotWidth,
        					plotHeight);
			}
	        break;
        case LAUNCH_SCATTER:
    		ext.setClipboard(sample+"\t"+loc);
    		ScatterPlot.createAndShowGUI(proj, new String[] {marker}, null, false);
	        break;
        case COPY_ID:
    		ext.setClipboard(sample);
	        break;
        case APPEND_ID:
    		ext.setClipboard(ext.getClipboard()+"\n"+sample);
	        break;
        case APPEND_ID_TO_FILE:
    		ext.setClipboard(sample);
    		ext.appendToFile(sample, filename);
	        break;
        case LABEL_ONLY:
	        break;
        default:
	        break;
        }
	}

	public boolean isEnabled() {
		switch (type) {
        case LAUNCH_TRAILER:
    		return Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY)+sample+Sample.SAMPLE_DATA_FILE_EXTENSION, jar); // needs to be updated anyway
        default:
        	return true;
        }
	}
}
