package hospelhornbg_svproject;

public interface OMIMUpdateListener {
	
	public void onDownloadStart();
	public void onDownloadProgressUpdate(int readBytes, int appTotalBytes);
	public void onDownloadFail();
	public void onDownloadComplete();
	
	public void onInvalidGeneSetFound();
	
	public void onTableReadStart();
	public void onReadTableLine(int lineNumber);
	public void onTableReadComplete();
	//public void onTableReadFail();
	
	public void onWritePrepareStart();
	public void onWritePrepareComplete(int geneCount);
	
	public void onWriteStart();
	public void onWriteGeneIDs(int geneNumber);
	public void onWriteComplete();
	
	
}
