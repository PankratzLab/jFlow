package org.genvisis.cnv.workflow;

import javax.swing.JComponent;

public interface StepAssist {

  /**
   * @return a JComponent that provides input assistance or suggestion for the implementing Step.
   *         This can be either a non-interactive component such as a JLabel, or could be a JButton
   *         that launches an entire subcomponent.
   */
  public JComponent getStepAssistComponent();

  public void updateStepAssistComponent(Variables variables);

}
