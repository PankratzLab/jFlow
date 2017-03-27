package org.genvisis.cnv.gui;

import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import net.miginfocom.swing.MigLayout;

public class JAccordionPanel extends JPanel {

	/**
	 * Auto-generated svUID
	 */
	private static final long serialVersionUID = 1L;

	private static final String UP = "-";
	private static final String DOWN = "+";

	private final JButton expandoButton;
	public JPanel topPanel;
	public JPanel contentPanel;
	private final JPanel panel;

	/**
	 * Create the panel.
	 */
	public JAccordionPanel() {
		setLayout(new MigLayout("hidemode 3", "[grow]", "[grow]"));

		topPanel = new JPanel();
		topPanel.setBorder(null);
		add(topPanel, "north");
		topPanel.setLayout(new MigLayout("", "[grow]", "[]"));

		panel = new JPanel();
		topPanel.add(panel, "east");
		panel.setBorder(null);
		panel.setLayout(new MigLayout("ins 4", "[]", "[]"));

		expandoButton = new JButton(UP);
		expandoButton.setFont(expandoButton.getFont().deriveFont(Font.BOLD, 14));
		panel.add(expandoButton, "cell 0 0");
		expandoButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (contentPanel.isVisible()) {
					shrink();
				} else {
					expand();
				}
			}
		});

		expandoButton.setMargin(new Insets(0, 5, 0, 5));
		expandoButton.setSelected(true);

		contentPanel = new JPanel();
		add(contentPanel, "cell 0 0,grow");
		contentPanel.setLayout(new MigLayout("", "[]", "[]"));
	}

	public void expand() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				contentPanel.setVisible(true);
				expandoButton.setText(UP);
				JAccordionPanel.this.invalidate();
				JAccordionPanel.this.repaint();
				for (JAccordionPanel jp : grp) {
					if (jp != JAccordionPanel.this) {
						jp.shrink();
					}
				}
			}
		});
	}

	public void shrink() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				contentPanel.setVisible(false);
				expandoButton.setText(DOWN);
				JAccordionPanel.this.invalidate();
				JAccordionPanel.this.repaint();
			}
		});
	}

	List<JAccordionPanel> grp = new ArrayList<JAccordionPanel>();

	public void addToGroup(List<JAccordionPanel> bg) {
		grp = bg;
		grp.add(this);
	}

}
