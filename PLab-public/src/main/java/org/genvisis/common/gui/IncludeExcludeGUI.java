package org.genvisis.common.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.HashMap;
import java.util.Map;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.pankratzlab.common.gui.UITools;

import edu.princeton.cs.algs4.RabinKarp;
import net.miginfocom.swing.MigLayout;

public class IncludeExcludeGUI extends JDialog {

  /**
  *
  */
  private static final long serialVersionUID = 1L;
  private final JPanel contentPanel = new JPanel();
  private volatile int close = JOptionPane.CLOSED_OPTION;
  private String[] opts;
  private Map<String, JCheckBox> optChks = new HashMap<>();
  private JTextField textField;
  private JList<JCheckBox> list;

  /**
   * Create the dialog.
   */
  public IncludeExcludeGUI(Frame owner, String[] opts, boolean[] preSel) {
    super(owner);
    initialize(opts, preSel);
  }

  /**
   * 
   * @param owner
   * @param opts
   * @param preSel
   * @wbp.parser.constructor
   */
  public IncludeExcludeGUI(Window owner, String[] opts, boolean[] preSel) {
    super(owner);
    initialize(opts, preSel);
  }

  private void initialize(String[] opts, boolean[] preSel) {
    setModal(true);
    setMinimumSize(new Dimension(100, 100));
    this.opts = opts;
    UITools.setSize(this, new Dimension(450, 300));
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    {
      DefaultListModel<JCheckBox> listModel = new DefaultListModel<>();
      for (int i = 0; i < opts.length; i++) {
        optChks.put(opts[i], new JCheckBox(opts[i], preSel[i]));
        listModel.addElement(optChks.get(opts[i]));
      }
      JScrollPane pane = new JScrollPane();
      list = new JList<>();
      list.setCellRenderer(new JCheckBoxListCellRenderer(getFont()));
      list.addMouseListener(new MouseAdapter() {

        int prevIndex = -1;

        public void mousePressed(MouseEvent e) {
          int index = list.locationToIndex(e.getPoint());
          if (index != -1) {
            JCheckBox checkbox = (JCheckBox) list.getModel().getElementAt(index);
            if (checkbox.isEnabled()) {
              boolean cl = e.getClickCount() >= 2;
              boolean li = prevIndex != -1 && index == prevIndex;
              if (cl || li) {
                checkbox.setSelected(!checkbox.isSelected());
              }
              prevIndex = index;
            }
            repaint();
          }
        }
      });
      list.addKeyListener(new KeyAdapter() {
        @Override
        public void keyPressed(KeyEvent e) {
          if (e.getKeyCode() == KeyEvent.VK_SPACE) {
            list.getSelectedValuesList().stream().forEach(chk -> {
              chk.setSelected(!chk.isSelected());
            });
          }
          super.keyPressed(e);
          repaint();
        }
      });
      contentPanel.setLayout(new MigLayout("", "[258px,grow]", "[17px][130px,grow][][]"));
      list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
      list.setModel(listModel);
      pane.setViewportView(list);
      contentPanel.add(pane, "cell 0 1,grow");
    }
    {
      JLabel lblSelectItemsTo = new JLabel("<html><u>Select Items to Include:</u></html>");
      lblSelectItemsTo.setFont(new Font("Arial", Font.BOLD, 14));
      contentPanel.add(lblSelectItemsTo, "cell 0 0,growx,aligny top");
    }
    {
      JLabel lblFind = new JLabel("Find:");
      contentPanel.add(lblFind, "cell 0 2");
    }
    {
      textField = new JTextField();
      textField.addKeyListener(
          new KeyAdapter() {

            @Override
            public void keyReleased(KeyEvent e) {
              super.keyTyped(e);
              SwingUtilities.invokeLater(
                  () -> {
                    String txt = textField.getText();
                    final DefaultListModel<JCheckBox> newMod = new DefaultListModel<>();
                    for (String opt : opts) {
                      if (new RabinKarp(txt).search(opt) < opt.length()) {
                        newMod.addElement(optChks.get(opt));
                      }
                    }
                    list.setModel(newMod);
                    list.revalidate();
                    repaint();
                  });
            }
          });
      contentPanel.add(textField, "cell 0 3,growx");
      textField.setColumns(10);
    }

    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        JButton okButton = new JButton("OK");
        okButton.addActionListener(new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent e) {
            close = JOptionPane.OK_OPTION;
            setVisible(false);
          }
        });
        okButton.setActionCommand("OK");
        buttonPane.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent e) {
            close = JOptionPane.CANCEL_OPTION;
            setVisible(false);
            dispose();
          }
        });
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
      }
    }
    pack();
  }

  /**
   * @return Either JOptionPane.OK_OPTION or JOptionPane.CANCEL_OPTION.
   */
  public int getCloseCode() {
    return close;
  }

  public boolean[] getSelected() {
    boolean[] sel = new boolean[opts.length];
    for (int i = 0; i < sel.length; i++) {
      sel[i] = optChks.get(opts[i]).isSelected();
    }
    return sel;
  }

}
