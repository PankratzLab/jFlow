package cnv.gui;

import cnv.plots.Trailer;
import common.ext;
import common.Files;
import java.awt.event.ActionEvent;
import javax.swing.AbstractAction;
import java.awt.Color;
import javax.swing.*;

import cnv.filesys.Project;

public class LaunchAction extends AbstractAction {
	public static final long serialVersionUID = 1L;
	public static final int LAUNCH_TRAILER = 1;
	public static final int COPY_ID = 2;
	public static final int APPEND_ID = 3;

	private Project proj;
	private String sample;
	private String filename;
	private String loc;
	private int type;
	private boolean jar;

	public LaunchAction(Project proj, String sample, String loc, Color color) {
		super(sample+" "+loc);
		this.type = LAUNCH_TRAILER;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		this.loc = loc;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}

	public LaunchAction(Project proj, String sample, Color color) {
		super(sample);
		this.type = COPY_ID;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}

	public LaunchAction(String filename, Project proj, String sample, Color color) {
		super(sample);
		this.type = APPEND_ID;
		this.filename = filename;
		this.proj = proj;
		this.jar = proj.getJarStatus();
		this.sample = sample;
		putValue(Action.SMALL_ICON, new ColorIcon(12, 12, color));
	}

	public void actionPerformed(ActionEvent e) {
		switch (type) {
        case LAUNCH_TRAILER:
    		ext.setClipboard(sample+"\t"+loc);
    		new Trailer(proj, sample, proj.getFilenames(Project.CNV_FILENAMES), loc.substring(0, loc.length()-1));
	        break;
        case COPY_ID:
    		ext.setClipboard(sample);
	        break;
        case APPEND_ID:
    		ext.setClipboard(sample);
    		ext.appendToFile(sample, filename);
	        break;
        default:
	        break;
        }
	}

	public boolean isEnabled() {
		switch (type) {
        case LAUNCH_TRAILER:
    		return Files.exists(proj.getDir(Project.IND_DIRECTORY)+sample+".samp", jar);
        default:
        	return true;
        }
	}
}
