package org.genvisis.cnv.filesys;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JButton;

import org.genvisis.cnv.Launch;
import org.genvisis.common.Grafik;

/**
 * Helper button configured to open a {@link ProjectPropertiesEditor} with a subset of properties.
 * Use this when you have a UI that depends on project properties that a user may want to change.
 */
public class PropertyEditorButton extends JButton {

	private static final long serialVersionUID = -2640225291773146965L;

	private final transient Project proj;
	private final String[] keys;

	/**
	 * Construct an editor button for the given project. The list of property keys will be applied as
	 * a filter in {@link ProjectPropertiesEditor#ProjectPropertiesEditor(Project, String...)}
	 */
	public PropertyEditorButton(Project project, String... projectPropertyKeys) {
		proj = project;
		keys = projectPropertyKeys;
		setIcon(Grafik.getImageIcon(ProjectPropertiesEditor.ICON));
		setPreferredSize(new Dimension(15, 15));
		setToolTipText("Edit related project properties");

		addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				//NB: copied the window logic from cnv.Launnch
				proj.getLog().report("Launching project properties editor...");
				final ProjectPropertiesEditor editor = new ProjectPropertiesEditor(proj, keys);
				editor.addWindowListener(new WindowAdapter() {
					@Override
					public void windowClosed(WindowEvent e) {
						Launch.getWindows()[0].requestFocus();
						editor.dispose();
					}
				});
				editor.setVisible(true);
			}
		});
	}


}
