package org.genvisis.fcs.sub;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.text.AttributedCharacterIterator;
import java.text.FieldPosition;
import java.text.Format;
import java.text.NumberFormat;
import java.text.ParseException;
import java.text.ParsePosition;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.genvisis.fcs.JFlow;
import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.Gate.RectangleGate;
import org.genvisis.fcs.gating.GateDimension.RectangleGateDimension;

import net.miginfocom.swing.MigLayout;

public class RectangleGateEditor extends JDialog {

  /** */
  private static final long serialVersionUID = 1L;

  private JFlow plot;
  private Gate gate;
  private final JPanel contentPanel = new JPanel();
  private JFormattedTextField txtXMin;
  private JFormattedTextField txtYMin;
  private JFormattedTextField txtXMax;
  private JFormattedTextField txtYMax;
  private JTextField txtName;
  private JLabel lblID;
  private JCheckBox chkXMaxUnbnd;
  private JCheckBox chkYMaxUnbnd;
  private JCheckBox chkYMinUnbnd;
  private JCheckBox chkXMinUnbnd;
  private JLabel lblYAxis;
  private JLabel lblXAxis;

  private volatile boolean cancelled = true;
  private final Format inputFormat = new ParseAllFormat(NumberFormat.getNumberInstance());

  /** Launch the application. */
  public static void main(String[] args) {
    try {
      RectangleGateEditor dialog = new RectangleGateEditor();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /** Create the dialog. */
  public RectangleGateEditor() {
    setTitle("Edit Gate");
    setModal(true);
    setModalityType(ModalityType.APPLICATION_MODAL);
    setBounds(100, 100, 450, 375);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(
        new MigLayout("", "[][grow][grow][grow][]", "[][][][][][][][][][][][][][]"));
    {
      JLabel lbl = new JLabel("ID:");
      contentPanel.add(lbl, "flowx,cell 1 0 3 1,alignx center");
    }
    {
      JLabel lbl = new JLabel("Name:");
      contentPanel.add(lbl, "flowx,cell 1 2 3 1,alignx center");
    }
    {
      JSeparator separator = new JSeparator();
      contentPanel.add(separator, "cell 1 4 3 1,growx");
    }
    {
      JLabel lbl = new JLabel("<html><u>X-Axis</u></html>");
      contentPanel.add(lbl, "cell 1 5,alignx center");
    }
    {
      JLabel lbl = new JLabel("<html><u>Y-Axis</u></html>");
      contentPanel.add(lbl, "cell 3 5,alignx center");
    }
    {
      lblXAxis = new JLabel("<AXIS PLACEHOLDER>");
      contentPanel.add(lblXAxis, "cell 1 6,alignx center");
    }
    {
      JSeparator separator = new JSeparator();
      separator.setOrientation(SwingConstants.VERTICAL);
      contentPanel.add(separator, "cell 2 5 1 8,alignx center,growy");
    }
    {
      lblYAxis = new JLabel("<AXIS PLACEHOLDER>");
      contentPanel.add(lblYAxis, "cell 3 6,alignx center");
    }
    {
      JLabel lbl = new JLabel("min:");
      contentPanel.add(lbl, "cell 1 7");
    }
    {
      JLabel lbl = new JLabel("min:");
      contentPanel.add(lbl, "cell 3 7");
    }
    {
      txtXMin = new JFormattedTextField(inputFormat);
      txtXMin.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
      txtXMin.addFocusListener(new MousePositionCorrectorListener());
      contentPanel.add(txtXMin, "cell 1 8,growx");
      txtXMin.setColumns(10);
    }
    {
      txtYMin = new JFormattedTextField(inputFormat);
      txtYMin.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
      txtYMin.addFocusListener(new MousePositionCorrectorListener());
      txtYMin.setColumns(10);
      contentPanel.add(txtYMin, "cell 3 8,growx");
    }
    {
      chkXMinUnbnd = new JCheckBox("unbounded");
      chkXMinUnbnd.addItemListener(
          new ItemListener() {

            @Override
            public void itemStateChanged(ItemEvent e) {
              txtXMin.setEnabled(e.getStateChange() == ItemEvent.DESELECTED);
            }
          });
      contentPanel.add(chkXMinUnbnd, "cell 1 9,alignx right");
    }
    {
      chkYMinUnbnd = new JCheckBox("unbounded");
      chkYMinUnbnd.addItemListener(
          new ItemListener() {

            @Override
            public void itemStateChanged(ItemEvent e) {
              txtYMin.setEnabled(e.getStateChange() == ItemEvent.DESELECTED);
            }
          });
      contentPanel.add(chkYMinUnbnd, "cell 3 9,alignx right");
    }
    {
      JLabel lbl = new JLabel("max:");
      contentPanel.add(lbl, "cell 1 10");
    }
    {
      JLabel lbl = new JLabel("max:");
      contentPanel.add(lbl, "cell 3 10");
    }
    {
      txtXMax = new JFormattedTextField(inputFormat);
      txtXMax.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
      txtXMax.addFocusListener(new MousePositionCorrectorListener());
      txtXMax.setColumns(10);
      contentPanel.add(txtXMax, "cell 1 11,growx");
    }
    {
      txtYMax = new JFormattedTextField(inputFormat);
      txtYMax.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
      txtYMax.addFocusListener(new MousePositionCorrectorListener());
      txtYMax.setColumns(10);
      contentPanel.add(txtYMax, "cell 3 11,growx");
    }
    {
      chkXMaxUnbnd = new JCheckBox("unbounded");
      chkXMaxUnbnd.addItemListener(
          new ItemListener() {

            @Override
            public void itemStateChanged(ItemEvent e) {
              txtXMax.setEnabled(e.getStateChange() == ItemEvent.DESELECTED);
            }
          });
      contentPanel.add(chkXMaxUnbnd, "cell 1 12,alignx right");
    }
    {
      chkYMaxUnbnd = new JCheckBox("unbounded");
      chkYMaxUnbnd.addItemListener(
          new ItemListener() {

            @Override
            public void itemStateChanged(ItemEvent e) {
              txtYMax.setEnabled(e.getStateChange() == ItemEvent.DESELECTED);
            }
          });
      contentPanel.add(chkYMaxUnbnd, "cell 3 12,alignx right");
    }
    {
      lblID = new JLabel("<AXIS_PLACEHOLDER>");
      contentPanel.add(lblID, "cell 1 0,alignx center");
    }
    {
      txtName = new JTextField();
      contentPanel.add(txtName, "cell 1 2,alignx center");
      txtName.setColumns(25);
    }
    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        JButton okButton = new JButton("OK");
        okButton.addActionListener(
            new ActionListener() {

              @Override
              public void actionPerformed(ActionEvent e) {
                String msg = valid();
                if (msg == null) {
                  cancelled = false;
                  setVisible(false);
                } else {
                  JOptionPane.showMessageDialog(
                      RectangleGateEditor.this, msg, "Error!", JOptionPane.ERROR_MESSAGE);
                }
              }
            });
        okButton.setActionCommand("OK");
        buttonPane.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(
            new ActionListener() {

              @Override
              public void actionPerformed(ActionEvent e) {
                setVisible(false);
                dispose();
              }
            });
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
      }
    }
  }

  public void setGate(JFlow plot, RectangleGate gate) {
    this.plot = plot;
    this.gate = gate;
    lblID.setText(gate.getID());
    txtName.setText(gate.getName());

    RectangleGateDimension xAxis = (RectangleGateDimension) gate.getXDimension();

    lblXAxis.setText(xAxis.getParam());

    boolean xMinUnbnd = !Float.isFinite(xAxis.getMin());
    chkXMinUnbnd.setSelected(xMinUnbnd);
    txtXMin.setEnabled(!xMinUnbnd);
    txtXMin.setText(xMinUnbnd ? "" : xAxis.getMin() + "");

    boolean xMaxUnbnd = !Float.isFinite(xAxis.getMax());
    chkXMaxUnbnd.setSelected(xMaxUnbnd);
    txtXMax.setEnabled(!xMaxUnbnd);
    txtXMax.setText(xMaxUnbnd ? "" : xAxis.getMax() + "");

    if (gate.getYDimension() != null) {
      RectangleGateDimension yAxis = (RectangleGateDimension) gate.getYDimension();
      lblYAxis.setText(yAxis.getParam());

      boolean yMinUnbnd = !Float.isFinite(yAxis.getMin());
      chkYMinUnbnd.setSelected(yMinUnbnd);
      chkYMinUnbnd.setEnabled(true);
      txtYMin.setEnabled(!yMinUnbnd);
      txtYMin.setText(yMinUnbnd ? "" : yAxis.getMin() + "");

      boolean yMaxUnbnd = !Float.isFinite(yAxis.getMax());
      chkYMaxUnbnd.setEnabled(true);
      chkYMaxUnbnd.setSelected(yMaxUnbnd);
      txtYMax.setEnabled(!yMaxUnbnd);
      txtYMax.setText(yMaxUnbnd ? "" : yAxis.getMax() + "");
    } else {
      lblYAxis.setText("");

      chkYMinUnbnd.setSelected(false);
      chkYMinUnbnd.setEnabled(false);
      txtYMin.setEnabled(false);
      txtYMin.setText("");

      chkYMaxUnbnd.setSelected(false);
      chkYMaxUnbnd.setEnabled(false);
      txtYMax.setEnabled(false);
      txtYMax.setText("");
    }
  }

  private String valid() {
    if ("".equals(getName())) {
      return "Name must not be blank.";
    }
    if (!gate.getName().equals(getName()) && plot.duplicateGateName(getName())) {
      return "New name must be unique.";
    }
    if (Float.isFinite(getXMin()) && Float.isFinite(getXMax()) && getXMin() > getXMax()) {
      return "X-min must be less than X-max.";
    }
    if (Float.isFinite(getYMin()) && Float.isFinite(getYMax()) && getYMin() > getYMax()) {
      return "Y-min must be less than Y-max.";
    }
    return null;
  }

  public boolean isCancelled() {
    return cancelled;
  }

  @Override
  public String getName() {
    return txtName.getText();
  }

  public float getXMin() {
    return chkXMinUnbnd.isSelected()
        ? Float.NEGATIVE_INFINITY
        : Float.parseFloat(txtXMin.getText());
  }

  public float getXMax() {
    return chkXMaxUnbnd.isSelected()
        ? Float.POSITIVE_INFINITY
        : Float.parseFloat(txtXMax.getText());
  }

  public float getYMin() {
    return chkYMinUnbnd.isSelected()
        ? Float.NEGATIVE_INFINITY
        : Float.parseFloat(txtYMin.getText());
  }

  public float getYMax() {
    return chkYMaxUnbnd.isSelected()
        ? Float.POSITIVE_INFINITY
        : Float.parseFloat(txtYMax.getText());
  }

  // https://stackoverflow.com/questions/1313390/is-there-any-way-to-accept-only-numeric-values-in-a-jtextfield
  // See http://tips4java.wordpress.com/2010/02/21/formatted-text-field-tips/
  private static class MousePositionCorrectorListener extends FocusAdapter {

    @Override
    public void focusGained(FocusEvent e) {
      /*
       * After a formatted text field gains focus, it replaces its text with its current value,
       * formatted appropriately of course. It does this after any focus listeners are notified. We
       * want to make sure that the caret is placed in the correct position rather than the default
       * that is before the 1st character!
       */
      final JTextField field = (JTextField) e.getSource();
      final int dot = field.getCaret().getDot();
      final int mark = field.getCaret().getMark();
      if (field.isEnabled() && field.isEditable()) {
        SwingUtilities.invokeLater(
            new Runnable() {

              @Override
              public void run() {
                // Only set the caret if the textfield hasn't got a selection on it
                if (dot == mark) {
                  field.getCaret().setDot(dot);
                }
              }
            });
      }
    }
  }

  /**
   * Decorator for a {@link Format Format} which only accepts values which can be completely parsed
   * by the delegate format. If the value can only be partially parsed, the decorator will refuse to
   * parse the value.
   */
  class ParseAllFormat extends Format {

    /** */
    private static final long serialVersionUID = 1L;

    private final Format fDelegate;

    /**
     * Decorate <code>aDelegate</code> to make sure if parser everything or nothing
     *
     * @param aDelegate The delegate format
     */
    public ParseAllFormat(Format aDelegate) {
      fDelegate = aDelegate;
    }

    @Override
    public StringBuffer format(Object obj, StringBuffer toAppendTo, FieldPosition pos) {
      return fDelegate.format(obj, toAppendTo, pos);
    }

    @Override
    public AttributedCharacterIterator formatToCharacterIterator(Object obj) {
      return fDelegate.formatToCharacterIterator(obj);
    }

    @Override
    public Object parseObject(String source, ParsePosition pos) {
      int initialIndex = pos.getIndex();
      Object result = fDelegate.parseObject(source, pos);
      if (result != null && pos.getIndex() < source.length()) {
        int errorIndex = pos.getIndex();
        pos.setIndex(initialIndex);
        pos.setErrorIndex(errorIndex);
        return null;
      }
      return result;
    }

    @Override
    public Object parseObject(String source) throws ParseException {
      // no need to delegate the call, super will call the parseObject( source, pos ) method
      return super.parseObject(source);
    }
  }
}
