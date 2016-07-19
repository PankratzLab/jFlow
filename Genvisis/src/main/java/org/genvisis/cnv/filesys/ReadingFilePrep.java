package cnv.filesys;

public interface ReadingFilePrep {
	public void init();

	public boolean validate();

	public void close();
}
