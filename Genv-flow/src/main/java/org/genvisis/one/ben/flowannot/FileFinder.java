package org.genvisis.one.ben.flowannot;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.math.BigInteger;
import java.util.List;
import java.util.Random;

import javax.swing.AbstractAction;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;
import scala.actors.threadpool.Arrays;

public class FileFinder extends JDialog {

	private final JPanel contentPanel = new JPanel();
	private JTextField textField;
	private String[] fullOptions;
	JList<String> list;

	public static List<String> showFileFinder(String[] options, boolean multiSelect) {
		FileFinder ff = new FileFinder(options, multiSelect);
		ff.setModal(true);
		ff.setVisible(true);
		return ff.list.getSelectedValuesList();
	}

	/**
	 * Create the dialog.
	 */
	private FileFinder(String[] options, boolean multiSelect) {
		this.fullOptions = options;
		setTitle("File Finder");
		setBounds(100, 100, 450, 300);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new MigLayout("", "[grow]", "[][][][grow]"));
		{
			JLabel lblSearch = new JLabel("Search:");
			contentPanel.add(lblSearch, "cell 0 0");
		}
		{
			textField = new JTextField();
			textField.addKeyListener(new KeyAdapter() {
				@Override
				public void keyReleased(KeyEvent e) {
					super.keyTyped(e);
					SwingUtilities.invokeLater(() -> {
						String txt = textField.getText();
						final DefaultListModel<String> newMod = new DefaultListModel<>();
						for (String opt : fullOptions) {
							if (RabinKarp.runRabinKarp(txt, opt) < opt.length()) {
								newMod.addElement(opt);
							}
						}
						list.setModel(newMod);
						list.revalidate();
						repaint();
					});
				}
			});
			contentPanel.add(textField, "cell 0 1,growx");
			textField.setColumns(10);
		}
		{
			JLabel lblFound = new JLabel("Found:");
			contentPanel.add(lblFound, "cell 0 2");
		}
		{
			JScrollPane scrollPane = new JScrollPane();
			contentPanel.add(scrollPane, "cell 0 3,grow");
			{
				list = new JList<String>((String[]) Arrays.copyOf(options, options.length));
				list.setSelectionMode(multiSelect ? ListSelectionModel.MULTIPLE_INTERVAL_SELECTION
																				 : ListSelectionModel.SINGLE_SELECTION);
				scrollPane.setViewportView(list);
			}
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				JButton okButton = new JButton();
				okButton.setAction(new AbstractAction() {
					@Override
					public void actionPerformed(ActionEvent e) {
						FileFinder.this.setVisible(false);
					}
				});
				okButton.setText("OK");
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton();
				cancelButton.setAction(new AbstractAction() {
					@Override
					public void actionPerformed(ActionEvent e) {
						FileFinder.this.list.clearSelection();
						FileFinder.this.setVisible(false);
					}
				});
				cancelButton.setText("Cancel");
				buttonPane.add(cancelButton);
			}
		}
	}


	/**
	 * Adapted into a static utility class (by cole0482) from
	 * https://algs4.cs.princeton.edu/53substring/RabinKarp.java.html
	 * 
	 * Copyright © 2000–2017, Robert Sedgewick and Kevin Wayne.
	 */
	public static class RabinKarp {

		public static int runRabinKarp(String pat, String txt) {
			int R = 256;
			int m = pat.length();
			long q = longRandomPrime();

			// precompute R^(m-1) % q for use in removing leading digit
			long RM = 1;
			for (int i = 1; i <= m - 1; i++)
				RM = (R * RM) % q;
			long patHash = hash(pat, m, R, q);

			int n = txt.length();
			if (n < m)
				return n;
			long txtHash = hash(txt, m, R, q);

			// check for match at offset 0
			if ((patHash == txtHash) && check(txt, 0, m, pat))
				return 0;

			// check for hash match; if hash match, check for exact match
			for (int i = m; i < n; i++) {
				// Remove leading digit, add trailing digit, check for match.
				txtHash = (txtHash + q - RM * txt.charAt(i - m) % q) % q;
				txtHash = (txtHash * R + txt.charAt(i)) % q;

				// match
				int offset = i - m + 1;
				if ((patHash == txtHash) && check(txt, offset, m, pat))
					return offset;
			}

			// no match
			return n;
		}

		// Compute hash for key[0..m-1].
		private static long hash(String key, int m, int R, long q) {
			long h = 0;
			for (int j = 0; j < m; j++)
				h = (R * h + key.charAt(j)) % q;
			return h;
		}

		private static boolean check(String txt, int i, int m, String pat) {
			for (int j = 0; j < m; j++)
				if (pat.charAt(j) != txt.charAt(i + j))
					return false;
			return true;
		}

		// a random 31-bit prime
		private static long longRandomPrime() {
			BigInteger prime = BigInteger.probablePrime(31, new Random());
			return prime.longValue();
		}
	}

}
