package org.seqcode.projects.seqview.model;

import java.io.IOException;

import org.seqcode.data.readdb.Client;
import org.seqcode.data.readdb.ClientException;

/**
 * All SeqViewModels that use ReadDB Clients should extend this class. 
 * Handles persistent connections correctly, including closing them after too long without comms. 
 * @author mahony
 *
 */
public abstract class ReadDBSeqViewModel extends SeqViewModel{

	protected Client client=null; 
	protected Long timeConnectionOpened;
	private Thread connectionTimeThread=null;
	private boolean connected=false;
	private final long CONNECTIONOPENLIMIT=5*60*1000; //Limit that connection can be open without comms = 5 minutes
	private final long CONNECTIONCHECKTIME=1*60*1000; //Check how long connection has been open every 1 minute
	
			
	public ReadDBSeqViewModel(){
		try {
			client = new Client(true);
			timeConnectionOpened = System.currentTimeMillis();
			connectionTimeThread = new Thread(new ReadDBSeqViewModelConnectionThread());
			connectionTimeThread.start();
		} catch (IOException | ClientException e) {
			e.printStackTrace();	
		}
	}
	
	public void openConnection(){
		synchronized(timeConnectionOpened){
			if(!connected)
				client.setPersistentConnection(true);
			connected=true;
			timeConnectionOpened = System.currentTimeMillis();
		}
	}
	
	public void close(){
		connectionTimeThread.interrupt();		
		if(client!=null){
			client.close();
			client=null;
		}
	}
	
	
	/**
	 * Simple thread to close persistent ReadDB connection after time limit of no comms
	 */
	class ReadDBSeqViewModelConnectionThread implements Runnable {	
		
		public void run(){
			while(true){
				try{
					Thread.sleep(CONNECTIONCHECKTIME);
					if(connected && System.currentTimeMillis() - timeConnectionOpened > CONNECTIONOPENLIMIT){
						synchronized(timeConnectionOpened){
							client.setPersistentConnection(false);
							connected=false;
						}
					}
				}catch(InterruptedException e){
					break;
				}
			}
		}
	}
}
