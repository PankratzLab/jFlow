package org.genvisis.flowannot;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Arrays;
import java.util.List;

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

import edu.princeton.cs.algs4.RabinKarp;
import net.miginfocom.swing.MigLayout;

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
              if (new RabinKarp(txt).search(opt) < opt.length()) {
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
        list = new JList<>((String[]) Arrays.copyOf(options, options.length));
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

}
