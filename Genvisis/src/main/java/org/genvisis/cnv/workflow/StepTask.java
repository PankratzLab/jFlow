package org.genvisis.cnv.workflow;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.pankratzlab.shared.gui.Task;

public class StepTask extends Task<Void, Void> {

  GenvisisWorkflowGUI gui;
  protected Throwable failureException;
  private ArrayList<String> failReasons = new ArrayList<>();
  protected Step.FINAL_CODE returnCode = Step.FINAL_CODE.CANCELLED;
  private boolean failed = false;
  private Step step;
  private List<Step> selectedSteps;
  private Variables variables;
  private Thread bgThread;

  public StepTask(GenvisisWorkflowGUI gui, Step step, List<Step> selectedSteps,
                  Variables variables) {
    this(gui, step, selectedSteps, variables, 0);
  }

  public StepTask(GenvisisWorkflowGUI gui, Step step, List<Step> selectedSteps, Variables variables,
                  int numUpdates) {
    super(step.getName(), numUpdates);
    this.gui = gui;
    this.step = step;
    this.selectedSteps = selectedSteps;
    this.variables = variables;
  }

  @Override
  protected Void doInBackground() throws Exception {
    this.gui.startStep(this.step);
    this.bgThread = Thread.currentThread();
    Exception e = null;
    Step.FINAL_CODE code = Step.FINAL_CODE.COMPLETE;
    try {
      this.step.setNecessaryPreRunProperties(variables);
      this.step.run(variables);
    } catch (RuntimeException e1) {
      if (e1.getCause() instanceof InterruptedException) {
        code = Step.FINAL_CODE.CANCELLED;
      }
      e = e1;
    } catch (Exception e1) {
      e = e1;
    }
    if (code != Step.FINAL_CODE.CANCELLED
        && (e != null || getFailed() || !this.step.checkIfOutputExists(variables))) {
      code = Step.FINAL_CODE.FAILED;
    }
    failureException = e;
    returnCode = code;
    return null;
  }

  @SuppressWarnings("deprecation")
  // uses Thread.stop();
  @Override
  protected void done() {
    if (this.isCancelled()) {
      bgThread.stop();
      returnCode = Step.FINAL_CODE.CANCELLED;
    }
    gui.endStep(this.step, returnCode);
    gui.nextStep(this, returnCode, selectedSteps, variables);
  }

  public void setFailed(String reason) {
    failed = true;
    failReasons.add(reason);
  }

  public List<String> getFailureMessages() {
    return failReasons;
  }

  public Throwable getFailureException() {
    return failureException;
  }

  public Step getStep() {
    return this.step;
  }

  public boolean getFailed() {
    return failed;
  }

}
