package cnv.manage;

import cnv.filesys.Project;

public interface TextExport {
    public void exportToText(Project proj, String outputFile);
}