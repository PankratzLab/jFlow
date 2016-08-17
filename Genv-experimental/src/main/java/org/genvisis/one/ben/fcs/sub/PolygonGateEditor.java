package org.genvisis.one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.text.AttributedCharacterIterator;
import java.text.FieldPosition;
import java.text.Format;
import java.text.NumberFormat;
import java.text.ParseException;
import java.text.ParsePosition;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.GateDimension;

import net.miginfocom.swing.MigLayout;

public class PolygonGateEditor extends JDialog {


  /**
   * <p>
   * Decorator for a {@link Format Format} which only accepts values which can be completely parsed
   * by the delegate format. If the value can only be partially parsed, the decorator will refuse to
   * parse the value.
   * </p>
   */
  class ParseAllFormat extends Format {
    /**
    * 
    */
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
    public Object parseObject(String source) throws ParseException {
      // no need to delegate the call, super will call the parseObject( source, pos ) method
      return super.parseObject(source);
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
  }

  /**
  * 
  */
  private static final long serialVersionUID = 1L;
  private final FCSPlot plot;
  private final Gate gate;
  private final JPanel contentPanel = new JPanel();
  private JTextField txtName;
  private JLabel lblID;
  private JLabel lblYAxis;

  private JLabel lblXAxis;
  private volatile boolean cancelled = true;
  private final Format inputFormat = new ParseAllFormat(NumberFormat.getNumberInstance());
  ArrayList<JFormattedTextField> xFields = new ArrayList<JFormattedTextField>();
  ArrayList<JFormattedTextField> yFields = new ArrayList<JFormattedTextField>();

  private JCheckBox chkMimicFlowJo;

  /**
   * Create the dialog.
   */
  public PolygonGateEditor(FCSPlot plot, PolygonGate gate) {
    this.plot = plot;
    this.gate = gate;
    setTitle("Edit Gate");
    setModal(true);
    setModalityType(ModalityType.APPLICATION_MODAL);
    setBounds(100, 100, 450, 375);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel
        .setLayout(new MigLayout("", "[][grow][][grow][]", "[][][][][][grow][][][][][][][]"));
    {
      JLabel lbl = new JLabel("ID:");
      contentPanel.add(lbl, "flowx,cell 1 0 3 1,alignx center");
    }
    {
      JLabel lbl = new JLabel("Name:");
      contentPanel.add(lbl, "flowx,cell 1 1 3 1,alignx center");
    }
    {
      JSeparator separator = new JSeparator();
      contentPanel.add(separator, "cell 1 2 3 1,growx");
    }
    {
      JLabel lbl = new JLabel("<html><u>X-Axis</u></html>");
      contentPanel.add(lbl, "cell 1 3,alignx center");
    }
    {
      JLabel lbl = new JLabel("<html><u>Y-Axis</u></html>");
      contentPanel.add(lbl, "cell 3 3,alignx center");
    }
    {
      GateDimension xAxis = gate.getDimensions().get(0);
      lblXAxis = new JLabel(xAxis.getParam());
      contentPanel.add(lblXAxis, "cell 1 4,alignx center");
    }
    {
      GateDimension yAxis = gate.getDimensions().get(1);
      lblYAxis = new JLabel(yAxis.getParam());
      contentPanel.add(lblYAxis, "cell 3 4,alignx center");
    }
    {
      lblID = new JLabel(gate.getID());
      contentPanel.add(lblID, "cell 1 0 3 1,alignx center");
    }
    {
      txtName = new JTextField(gate.getName());
      contentPanel.add(txtName, "cell 1 1 3 1,alignx center");
      txtName.setColumns(25);
    }
    {
      JScrollPane scrollPane = new JScrollPane();
      contentPanel.add(scrollPane, "cell 1 5 3 7,grow");
      {
        JPanel panel = new JPanel();
        scrollPane.setViewportView(panel);
        panel.setLayout(new MigLayout("", "[grow][][grow]", "[][][]"));


        Path2D path = gate.getPath();
        PathIterator pi = path.getPathIterator(null);
        double[] coords = new double[6];
        int index = 0;
        while (!pi.isDone()) {
          int cd = pi.currentSegment(coords);
          if (cd == PathIterator.SEG_CLOSE) {
            pi.next();
            continue;
          }
          JFormattedTextField xField = new JFormattedTextField(inputFormat);
          xField.setText(coords[0] + "");
          JFormattedTextField yField = new JFormattedTextField(inputFormat);
          yField.setText(coords[1] + "");

          panel.add(xField, "cell 0 " + index + ", growx");
          panel.add(yField, "cell 2 " + index + ", growx");
          index++;
          xFields.add(xField);
          yFields.add(yField);

          pi.next();
        }

        JSeparator separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator, "cell 1 0 1 " + index + ",growy");
      }
    }
    {
      chkMimicFlowJo = new JCheckBox("Mimic FlowJo Binning", gate.getMimicsFlowJo());
      contentPanel.add(chkMimicFlowJo, "cell 1 12 3 1,alignx right");
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
            String msg = valid();
            if (msg == null) {
              cancelled = false;
              setVisible(false);
            } else {
              JOptionPane.showMessageDialog(PolygonGateEditor.this, msg, "Error!",
                  JOptionPane.ERROR_MESSAGE);
            }
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
            setVisible(false);
            dispose();
          }
        });
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
      }
    }
  }

  public boolean getMimicFlowJo() {
    return chkMimicFlowJo.isSelected();
  }

  @Override
  public String getName() {
    return txtName.getText();
  }

  public Path2D getNewPath() {
    Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
    float x, y;
    x = Float.parseFloat(xFields.get(0).getValue().toString());
    y = Float.parseFloat(yFields.get(0).getValue().toString());
    path.moveTo(x, y);
    for (int i = 1; i < xFields.size(); ++i) {
      x = Float.parseFloat(xFields.get(i).getValue().toString());
      y = Float.parseFloat(yFields.get(i).getValue().toString());
      path.lineTo(x, y);
    }
    path.closePath();
    return path;
  }

  public boolean isCancelled() {
    return cancelled;
  }

  private String valid() {
    if ("".equals(getName())) {
      return "Name must not be blank.";
    }
    if (!gate.getName().equals(getName()) && plot.duplicateGateName(getName())) {
      return "New name must be unique.";
    }
    return null;
  }

}
